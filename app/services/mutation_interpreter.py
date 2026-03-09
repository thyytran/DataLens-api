# backend/services/mutation_interpreter.py

from Bio.PDB import PDBParser, DSSP
from Bio.PDB.SASA import ShrakeRupley
from typing import Dict, List, Optional
from pathlib import Path
from dataclasses import dataclass
from .amino_acid_properties import AMINO_ACIDS, compare_amino_acids


@dataclass
class FoldXResult:
    """Parsed FoldX output"""
    ddg: float
    sd: float
    wt_energy: float
    mut_energy: float
    num_runs: int
    energy_breakdown: Dict[str, float]
    
    @classmethod
    def from_raw_file(cls, raw_file: str):
        """Parse FoldX Raw output file"""
        with open(raw_file, 'r') as f:
            lines = f.readlines()
        
        # Find data lines
        data_lines = []
        for line in lines:
            if line.strip() and not line.startswith('FoldX') and \
               not line.startswith('by') and not line.startswith('Pdb') and \
               not line.startswith('---') and not line.startswith('PDB') and \
               not line.startswith('Output'):
                parts = line.split('\t')
                if len(parts) > 2:
                    data_lines.append(parts)
        
        # Separate mutant and wild-type
        mutant_lines = [l for l in data_lines if not l[0].startswith('WT_')]
        wt_lines = [l for l in data_lines if l[0].startswith('WT_')]
        
        if not mutant_lines or not wt_lines:
            raise ValueError("Could not find mutant or wild-type data in Raw file")
        
        # Average energies
        mut_energies = [float(l[1]) for l in mutant_lines]
        wt_energies = [float(l[1]) for l in wt_lines]
        
        avg_mut = sum(mut_energies) / len(mut_energies)
        avg_wt = sum(wt_energies) / len(wt_energies)
        
        ddg = avg_mut - avg_wt
        
        # Standard deviation
        import statistics
        sd_mut = statistics.stdev(mut_energies) if len(mut_energies) > 1 else 0.0
        sd_wt = statistics.stdev(wt_energies) if len(wt_energies) > 1 else 0.0
        combined_sd = (sd_mut**2 + sd_wt**2) ** 0.5
        
        # Energy breakdown (difference between first mutant and WT)
        energy_breakdown = {
            'backbone_hbond': float(mutant_lines[0][2]) - float(wt_lines[0][2]),
            'sidechain_hbond': float(mutant_lines[0][3]) - float(wt_lines[0][3]),
            'van_der_waals': float(mutant_lines[0][4]) - float(wt_lines[0][4]),
            'electrostatics': float(mutant_lines[0][5]) - float(wt_lines[0][5]),
            'solvation_polar': float(mutant_lines[0][6]) - float(wt_lines[0][6]),
            'solvation_hydrophobic': float(mutant_lines[0][7]) - float(wt_lines[0][7]),
            'clash': float(mutant_lines[0][8]) - float(wt_lines[0][8]),
            'entropy_sidechain': float(mutant_lines[0][9]) - float(wt_lines[0][9]),
            'entropy_mainchain': float(mutant_lines[0][10]) - float(wt_lines[0][10]),
            'torsional_clash': float(mutant_lines[0][13]) - float(wt_lines[0][13]),
            'backbone_clash': float(mutant_lines[0][14]) - float(wt_lines[0][14]),
        }
        
        return cls(
            ddg=ddg,
            sd=combined_sd,
            wt_energy=avg_wt,
            mut_energy=avg_mut,
            num_runs=len(mutant_lines),
            energy_breakdown=energy_breakdown
        )


class MutationInterpreter:
    """
    HOPE-style mutation interpretation system with FoldX integration
    Analyzes structural and biochemical impacts of mutations
    """
    
    def __init__(self, pdb_file: str, foldx_output_dir: Optional[str] = None):
        """
        Initialize with PDB structure file
        
        Args:
            pdb_file: Path to PDB file (preferably repaired)
            foldx_output_dir: Optional directory containing FoldX .fxout files
        """
        self.pdb_file = Path(pdb_file)
        self.foldx_dir = Path(foldx_output_dir) if foldx_output_dir else None
        
        self.parser = PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('protein', pdb_file)
        
        # Calculate SASA for all residues
        self.sasa_dict = self._calculate_sasa()
        
        # Get secondary structure assignments
        self.ss_dict = self._get_secondary_structure(pdb_file)
    
    
    def _calculate_sasa(self) -> Dict:
        """
        Calculate Solvent Accessible Surface Area for all residues
        Uses Shrake-Rupley algorithm via BioPython
        
        Returns:
            dict: {(chain, resnum): sasa_value}
        """
        sr = ShrakeRupley()
        sr.compute(self.structure, level="R")  # Residue level
        
        sasa_dict = {}
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] == ' ':  # Skip heteroatoms
                        key = (chain.id, residue.id[1])
                        sasa_dict[key] = residue.sasa
        
        return sasa_dict
    
    
    def _get_secondary_structure(self, pdb_file: str) -> Dict:
        """
        Get secondary structure assignments using DSSP
        
        Returns:
            dict: {(chain, resnum): ss_code}
        """
        try:
            model = self.structure[0]
            dssp = DSSP(model, pdb_file, dssp='mkdssp')
            
            ss_dict = {}
            for key in dssp:
                chain_id = key[0]
                res_id = key[1][1]
                ss_code = dssp[key][2]
                ss_dict[(chain_id, res_id)] = ss_code
            
            return ss_dict
            
        except Exception as e:
            print(f"Warning: Could not calculate secondary structure: {e}")
            return {}
    
    
    def load_foldx_results(self, mutation: str) -> Optional[FoldXResult]:
        """
        Load FoldX results for a mutation
        
        Args:
            mutation: Mutation string like "SA729A"
            
        Returns:
            FoldXResult or None if not found
        """
        if not self.foldx_dir:
            return None
        
        # Look for Raw file
        pdb_stem = self.pdb_file.stem
        raw_file = self.foldx_dir / f"Raw_{pdb_stem}.fxout"
        
        if not raw_file.exists():
            # Try alternative naming
            raw_file = self.foldx_dir / f"Raw_{pdb_stem}_Repair.fxout"
        
        if not raw_file.exists():
            print(f"Warning: FoldX Raw file not found at {raw_file}")
            return None
        
        try:
            return FoldXResult.from_raw_file(str(raw_file))
        except Exception as e:
            print(f"Error parsing FoldX results: {e}")
            return None
    
    
    def classify_location(self, chain: str, position: int) -> Dict:
        """
        Classify residue as core, surface, or partially buried
        
        Args:
            chain: Chain ID
            position: Residue number
            
        Returns:
            dict with location info
        """
        key = (chain, position)
        
        if key not in self.sasa_dict:
            return {
                'location': 'unknown',
                'sasa': None,
                'description': 'Unable to calculate SASA'
            }
        
        sasa = self.sasa_dict[key]
        
        # Classification thresholds (Å²)
        if sasa < 20:
            location = 'core'
            description = 'buried in protein core'
        elif sasa > 60:
            location = 'surface'
            description = 'exposed on protein surface'
        else:
            location = 'partially_buried'
            description = 'partially buried'
        
        return {
            'location': location,
            'sasa': sasa,
            'description': description
        }
    
    
    def get_secondary_structure(self, chain: str, position: int) -> Dict:
        """
        Get secondary structure assignment for residue
        
        Returns:
            dict with SS info
        """
        key = (chain, position)
        
        if key not in self.ss_dict:
            return {
                'ss_type': 'unknown',
                'description': 'Secondary structure not available'
            }
        
        ss_code = self.ss_dict[key]
        
        # Map DSSP codes to descriptions
        ss_map = {
            'H': {'type': 'helix', 'description': 'alpha helix'},
            'G': {'type': 'helix', 'description': '3-10 helix'},
            'I': {'type': 'helix', 'description': 'pi helix'},
            'E': {'type': 'sheet', 'description': 'beta sheet'},
            'B': {'type': 'sheet', 'description': 'beta bridge'},
            'T': {'type': 'turn', 'description': 'turn'},
            'S': {'type': 'bend', 'description': 'bend'},
            'C': {'type': 'coil', 'description': 'random coil'},
            '-': {'type': 'coil', 'description': 'random coil'}
        }
        
        ss_info = ss_map.get(ss_code, {'type': 'unknown', 'description': 'unassigned'})
        ss_info['code'] = ss_code
        
        return ss_info
    
    
    def analyze_mutation(
        self, 
        wt_residue: str, 
        mut_residue: str, 
        chain: str,
        position: int,
        alphamissense_score: Optional[float] = None,
        foldx_ddg: Optional[float] = None,
        foldx_results: Optional[FoldXResult] = None
    ) -> Dict:
        """
        Complete mutation analysis with HOPE-style interpretations
        
        Args:
            wt_residue: Wild-type amino acid (single letter)
            mut_residue: Mutant amino acid (single letter)
            chain: Chain ID
            position: Residue position
            alphamissense_score: Optional AlphaMissense pathogenicity score
            foldx_ddg: Optional FoldX ΔΔG (if already calculated)
            foldx_results: Optional parsed FoldX results object
            
        Returns:
            dict with comprehensive analysis
        """
        # Construct mutation string
        mutation = f"{wt_residue}{chain}{position}{mut_residue}"
        
        # Load FoldX results if available
        if foldx_results is None and self.foldx_dir:
            foldx_results = self.load_foldx_results(mutation)
        
        # Get amino acid properties comparison
        aa_comparison = compare_amino_acids(wt_residue, mut_residue)
        
        # Get structural context
        location_info = self.classify_location(chain, position)
        ss_info = self.get_secondary_structure(chain, position)
        
        # Generate interpretations
        interpretations = []
        
        # 1. CHARGE ANALYSIS
        interpretations.extend(
            self._analyze_charge(aa_comparison, location_info)
        )
        
        # 2. SIZE ANALYSIS
        interpretations.extend(
            self._analyze_size(aa_comparison, location_info, ss_info)
        )
        
        # 3. HYDROPHOBICITY ANALYSIS
        interpretations.extend(
            self._analyze_hydrophobicity(aa_comparison, location_info)
        )
        
        # 4. SECONDARY STRUCTURE SPECIFIC ANALYSIS
        interpretations.extend(
            self._analyze_secondary_structure(aa_comparison, ss_info)
        )
        
        # 5. SPECIAL RESIDUE ANALYSIS (Gly, Pro, Cys)
        interpretations.extend(
            self._analyze_special_residues(wt_residue, mut_residue, ss_info)
        )
        
        # 6. FOLDX ENERGY ANALYSIS (if available)
        if foldx_results:
            interpretations.extend(
                self._analyze_foldx_energy(foldx_results, location_info)
            )
            foldx_ddg = foldx_results.ddg
        
        # Generate summary
        summary = self._generate_summary(
            interpretations, 
            alphamissense_score, 
            foldx_ddg,
            foldx_results
        )
        
        result = {
            'mutation': mutation,
            'location': location_info,
            'secondary_structure': ss_info,
            'property_changes': aa_comparison['differences'],
            'interpretations': interpretations,
            'summary': summary,
            'predictions': {
                'alphamissense': alphamissense_score,
                'foldx_ddg': foldx_ddg
            }
        }
        
        # Add FoldX details if available
        if foldx_results:
            result['foldx_details'] = {
                'ddg': foldx_results.ddg,
                'sd': foldx_results.sd,
                'wt_energy': foldx_results.wt_energy,
                'mut_energy': foldx_results.mut_energy,
                'num_runs': foldx_results.num_runs,
                'energy_breakdown': foldx_results.energy_breakdown
            }
        
        return result
    
    
    def _analyze_charge(self, aa_comp: Dict, location: Dict) -> List[Dict]:
        """Analyze charge-related impacts"""
        interpretations = []
        charge_change = aa_comp['differences']['charge']
        
        if charge_change['changed']:
            if charge_change['type'] == 'charge_reversal':
                interpretations.append({
                    'category': 'charge',
                    'severity': 'high',
                    'icon': '⚠️',
                    'text': (
                        f"The mutation introduces opposite charge "
                        f"({aa_comp['reference']['charge']} → {aa_comp['alternative']['charge']}). "
                        f"This likely disrupts electrostatic interactions."
                    )
                })
                
                if location['location'] == 'surface':
                    interpretations.append({
                        'category': 'charge',
                        'severity': 'high',
                        'icon': '🔌',
                        'text': (
                            "Charge reversal on protein surface may disrupt "
                            "interactions with binding partners or DNA/RNA."
                        )
                    })
            
            elif charge_change['type'] == 'gain_charge':
                severity = 'medium' if location['location'] == 'surface' else 'high'
                interpretations.append({
                    'category': 'charge',
                    'severity': severity,
                    'icon': '➕',
                    'text': (
                        f"Introduction of charged residue "
                        f"({aa_comp['alternative']['charge']}) "
                        f"at {location['description']} position."
                    )
                })
            
            elif charge_change['type'] == 'loss_charge':
                interpretations.append({
                    'category': 'charge',
                    'severity': 'medium',
                    'icon': '➖',
                    'text': (
                        f"Loss of {aa_comp['reference']['charge']} charge "
                        "may affect local electrostatic environment."
                    )
                })
        
        return interpretations
    
    
    def _analyze_size(self, aa_comp: Dict, location: Dict, ss: Dict) -> List[Dict]:
        """Analyze size-related impacts"""
        interpretations = []
        size_diff = aa_comp['differences']['size']
        
        if abs(size_diff['value']) >= 2:
            if size_diff['description'] == 'larger':
                severity = 'high' if location['location'] == 'core' else 'medium'
                interpretations.append({
                    'category': 'size',
                    'severity': severity,
                    'icon': '📏',
                    'text': (
                        f"Larger residue ({aa_comp['alternative']['size']} vs "
                        f"{aa_comp['reference']['size']}) "
                        f"{'in protein core ' if location['location'] == 'core' else ''}"
                        "may cause steric clashes with neighboring residues."
                    )
                })
            
            elif size_diff['description'] == 'smaller':
                if location['location'] == 'core':
                    interpretations.append({
                        'category': 'size',
                        'severity': 'high',
                        'icon': '🕳️',
                        'text': (
                            f"Smaller residue in protein core "
                            f"({aa_comp['reference']['size']} → {aa_comp['alternative']['size']}) "
                            "may create destabilizing cavity, disrupting hydrophobic packing."
                        )
                    })
                elif location['location'] == 'surface':
                    interpretations.append({
                        'category': 'size',
                        'severity': 'low',
                        'icon': 'ℹ️',
                        'text': (
                            f"Smaller residue on surface may reduce "
                            "interactions with binding partners."
                        )
                    })
        
        return interpretations
    
    
    def _analyze_hydrophobicity(self, aa_comp: Dict, location: Dict) -> List[Dict]:
        """Analyze hydrophobicity changes"""
        interpretations = []
        hydro_change = aa_comp['differences']['hydrophobicity']
        
        if abs(hydro_change['value']) > 2.0:
            
            if (hydro_change['description'] == 'more_hydrophilic' and 
                location['location'] == 'core'):
                interpretations.append({
                    'category': 'hydrophobicity',
                    'severity': 'high',
                    'icon': '💧',
                    'text': (
                        "Introduction of hydrophilic residue in hydrophobic core "
                        "is highly destabilizing. This disrupts protein folding."
                    )
                })
            
            elif (hydro_change['description'] == 'more_hydrophobic' and 
                  location['location'] == 'surface'):
                interpretations.append({
                    'category': 'hydrophobicity',
                    'severity': 'medium',
                    'icon': '🛡️',
                    'text': (
                        "Introduction of hydrophobic residue on protein surface "
                        "may promote aggregation or alter solubility."
                    )
                })
            
            elif (hydro_change['description'] == 'more_hydrophilic' and 
                  location['location'] == 'surface'):
                interpretations.append({
                    'category': 'hydrophobicity',
                    'severity': 'low',
                    'icon': 'ℹ️',
                    'text': (
                        "Change to more hydrophilic residue on surface "
                        "is generally well-tolerated."
                    )
                })
        
        return interpretations
    
    
    def _analyze_secondary_structure(self, aa_comp: Dict, ss: Dict) -> List[Dict]:
        """Analyze secondary structure-specific impacts"""
        interpretations = []
        
        if ss['ss_type'] == 'unknown':
            return interpretations
        
        wt_name = aa_comp['reference']['name']
        mut_name = aa_comp['alternative']['name']
        
        if ss['ss_type'] == 'helix':
            if mut_name == 'Proline':
                interpretations.append({
                    'category': 'secondary_structure',
                    'severity': 'high',
                    'icon': '🔀',
                    'text': (
                        "Introduction of Proline into helix is severely "
                        "disruptive. Proline cannot adopt helical geometry."
                    )
                })
            
            elif mut_name == 'Glycine' and wt_name != 'Glycine':
                interpretations.append({
                    'category': 'secondary_structure',
                    'severity': 'medium',
                    'icon': '🌀',
                    'text': (
                        "Introduction of Glycine into helix increases "
                        "backbone flexibility and may destabilize structure."
                    )
                })
        
        elif ss['ss_type'] == 'sheet':
            if (aa_comp['reference']['size_value'] >= 3 and 
                aa_comp['alternative']['size_value'] <= 2):
                interpretations.append({
                    'category': 'secondary_structure',
                    'severity': 'medium',
                    'icon': '📋',
                    'text': (
                        f"Mutation to smaller residue in beta sheet may "
                        "affect strand packing and stability."
                    )
                })
        
        elif ss['ss_type'] in ['turn', 'coil']:
            if wt_name in ['Glycine', 'Proline'] and mut_name not in ['Glycine', 'Proline']:
                interpretations.append({
                    'category': 'secondary_structure',
                    'severity': 'medium',
                    'icon': '🔄',
                    'text': (
                        f"Loss of {wt_name} from flexible loop region may "
                        "reduce conformational flexibility."
                    )
                })
        
        return interpretations
    
    
    def _analyze_special_residues(self, wt: str, mut: str, ss: Dict) -> List[Dict]:
        """Analyze special cases: Glycine, Proline, Cysteine"""
        interpretations = []
        
        wt_props = AMINO_ACIDS[wt]
        mut_props = AMINO_ACIDS[mut]
        
        # Cysteine changes
        if wt_props['sulphur'] and not mut_props['sulphur']:
            interpretations.append({
                'category': 'special',
                'severity': 'high',
                'icon': '🔗',
                'text': (
                    "Loss of Cysteine may disrupt disulfide bond formation, "
                    "which is critical for protein stability."
                )
            })
        elif not wt_props['sulphur'] and mut_props['sulphur']:
            interpretations.append({
                'category': 'special',
                'severity': 'medium',
                'icon': '🔗',
                'text': (
                    "Introduction of Cysteine may form aberrant disulfide bonds "
                    "or interfere with existing bonds."
                )
            })
        
        # Glycine is special
        if wt == 'G' and mut != 'G':
            interpretations.append({
                'category': 'special',
                'severity': 'medium',
                'icon': '🎯',
                'text': (
                    "Loss of Glycine removes unique backbone flexibility. "
                    "Glycine is often critical in tight turns and compact structures."
                )
            })
        
        # Proline is special
        if wt == 'P' and mut != 'P':
            interpretations.append({
                'category': 'special',
                'severity': 'medium',
                'icon': '🔒',
                'text': (
                    "Loss of Proline increases backbone flexibility at this position."
                )
            })
        elif mut == 'P' and wt != 'P':
            interpretations.append({
                'category': 'special',
                'severity': 'medium',
                'icon': '🔒',
                'text': (
                    "Introduction of Proline restricts backbone conformation "
                    "and may disrupt existing structure."
                )
            })
        
        return interpretations
    
    
    def _analyze_foldx_energy(self, foldx: FoldXResult, location: Dict) -> List[Dict]:
        """
        Analyze FoldX energy breakdown to generate specific interpretations
        """
        interpretations = []
        breakdown = foldx.energy_breakdown
        
        # Significant threshold
        THRESHOLD = 0.1
        
        # Backbone H-bonds
        if breakdown['backbone_hbond'] > THRESHOLD:
            interpretations.append({
                'category': 'foldx_energy',
                'severity': 'medium',
                'icon': '🔗',
                'text': f"Loss of backbone hydrogen bonds (+{breakdown['backbone_hbond']:.2f} kcal/mol) reduces stability."
            })
        
        # Sidechain H-bonds
        if breakdown['sidechain_hbond'] > THRESHOLD:
            interpretations.append({
                'category': 'foldx_energy',
                'severity': 'high',
                'icon': '🔗',
                'text': f"Loss of sidechain hydrogen bonds (+{breakdown['sidechain_hbond']:.2f} kcal/mol) is destabilizing."
            })
        
        # Van der Waals
        if breakdown['van_der_waals'] > THRESHOLD:
            interpretations.append({
                'category': 'foldx_energy',
                'severity': 'medium',
                'icon': '⚛️',
                'text': f"Loss of van der Waals contacts (+{breakdown['van_der_waals']:.2f} kcal/mol) weakens structure."
            })
        
        # Clashes
        if breakdown['clash'] > THRESHOLD:
            interpretations.append({
                'category': 'foldx_energy',
                'severity': 'high',
                'icon': '💥',
                'text': f"Steric clashes introduced (+{breakdown['clash']:.2f} kcal/mol). Mutant residue doesn't fit well."
            })
        
        # Hydrophobic solvation
        if breakdown['solvation_hydrophobic'] > THRESHOLD:
            interpretations.append({
                'category': 'foldx_energy',
                'severity': 'medium',
                'icon': '💧',
                'text': f"Poor hydrophobic burial (+{breakdown['solvation_hydrophobic']:.2f} kcal/mol) increases energy."
            })
        
        # Entropy changes
        total_entropy = breakdown['entropy_sidechain'] + breakdown['entropy_mainchain']
        if total_entropy > THRESHOLD:
            interpretations.append({
                'category': 'foldx_energy',
                'severity': 'low',
                'icon': '♻️',
                'text': f"Loss of conformational entropy (+{total_entropy:.2f} kcal/mol) contributes to destabilization."
            })
        
        return interpretations
    
    
    def _generate_summary(
        self, 
        interpretations: List[Dict],
        alphamissense: Optional[float],
        foldx_ddg: Optional[float],
        foldx_results: Optional[FoldXResult]
    ) -> str:
        """Generate HOPE-style summary text"""
        
        # Count severity levels
        high_severity = sum(1 for i in interpretations if i['severity'] == 'high')
        
        # Build summary
        summary_parts = []
        
        # Overall assessment
        if high_severity >= 2:
            summary_parts.append(
                "This mutation is predicted to be highly damaging to protein structure and function."
            )
        elif high_severity == 1:
            summary_parts.append(
                "This mutation is likely to impact protein structure or function."
            )
        else:
            summary_parts.append(
                "This mutation may have mild to moderate effects on protein structure."
            )
        
        # Add prediction scores
        if alphamissense is not None:
            if alphamissense > 0.7:
                summary_parts.append(
                    f"AlphaMissense predicts this as pathogenic (score: {alphamissense:.3f})."
                )
            elif alphamissense < 0.3:
                summary_parts.append(
                    f"AlphaMissense predicts this as benign (score: {alphamissense:.3f})."
                )
            else:
                summary_parts.append(
                    f"AlphaMissense gives an ambiguous prediction (score: {alphamissense:.3f})."
                )
        
        if foldx_ddg is not None:
            if foldx_ddg > 2.0:
                summary_parts.append(
                    f"FoldX predicts significant destabilization (ΔΔG: +{foldx_ddg:.2f} kcal/mol)."
                )
            elif foldx_ddg > 1.0:
                summary_parts.append(
                    f"FoldX predicts moderate destabilization (ΔΔG: +{foldx_ddg:.2f} kcal/mol)."
                )
            elif foldx_ddg < -1.0:
                summary_parts.append(
                    f"FoldX predicts stabilization (ΔΔG: {foldx_ddg:.2f} kcal/mol)."
                )
            else:
                summary_parts.append(
                    f"FoldX predicts minimal stability change (ΔΔG: {foldx_ddg:.2f} kcal/mol)."
                )
            
            # Add energy insights if available
            if foldx_results:
                breakdown = foldx_results.energy_breakdown
                
                # Find main contributor
                main_term = max(breakdown.items(), key=lambda x: abs(x[1]))
                if abs(main_term[1]) > 0.1:
                    term_name = main_term[0].replace('_', ' ')
                    if main_term[1] > 0:
                        summary_parts.append(
                            f"The destabilization is primarily driven by {term_name} "
                            f"(+{main_term[1]:.2f} kcal/mol)."
                        )
        
        return " ".join(summary_parts)


# ==================== Usage Example ====================
"""
if __name__ == "__main__":
    # Initialize with PDB and FoldX output directory
    interpreter = MutationInterpreter(
        pdb_file="/path/to/7lmk_Repair.pdb",
        foldx_output_dir="/path/to/foldx_output"
    )
    
    # Analyze mutation SA729A
    result = interpreter.analyze_mutation(
        wt_residue='S',
        mut_residue='A',
        chain='A',
        position=729,
        alphamissense_score=0.12
        # foldx_ddg will be loaded automatically from output files
    )
    
    # Print results
    print(f"\n{result['mutation']}")
    print(f"Location: {result['location']['description']} (SASA: {result['location']['sasa']:.1f} Ų)")
    print(f"Secondary structure: {result['secondary_structure']['description']}")
    
    if 'foldx_details' in result:
        foldx = result['foldx_details']
        print(f"\nFoldX Analysis:")
        print(f"  ΔΔG: {foldx['ddg']:+.2f} ± {foldx['sd']:.2f} kcal/mol ({foldx['num_runs']} runs)")
        print(f"  Top energy contributors:")
        for term, value in sorted(foldx['energy_breakdown'].items(), 
                                   key=lambda x: abs(x[1]), reverse=True)[:3]:
            if abs(value) > 0.05:
                print(f"    {term}: {value:+.3f} kcal/mol")
    
    print(f"\nSummary: {result['summary']}")
    print("\nDetailed interpretations:")
    for interp in result['interpretations']:
        print(f"  {interp['icon']} [{interp['severity'].upper()}] {interp['text']}")
        
        
"""