# backend/services/datalens_report.py

from typing import Dict, List, Optional
from dataclasses import dataclass
from enum import Enum

# Import your other modules
from .amino_acid_properties import AMINO_ACIDS, compare_amino_acids
from .mutation_interpreter import MutationInterpreter


class Severity(Enum):
    HIGH = "high"
    MEDIUM = "medium"
    LOW = "low"
    BENIGN = "benign"


@dataclass
class Interpretation:
    """Single interpretation finding"""
    category: str  # 'charge', 'size', 'hydrophobicity', 'secondary_structure', etc.
    severity: Severity
    icon: str
    text: str
    evidence: Optional[Dict] = None


class DataLensReport:
    """
    Complete mutation analysis report combining:
    - FoldX stability predictions (ΔΔG)
    - AlphaMissense pathogenicity scores
    - Structural context (SASA, SS, location)
    - Amino acid property analysis
    - HOPE-style interpretations
    """
    
    def __init__(
        self,
        mutation: str,  # e.g., "SA729A"
        pdb_id: str,
        chain: str,
        position: int,
        wt_aa: str,
        mut_aa: str,
        foldx_results: Dict,
        alphamissense_score: Optional[float] = None,
        structure_context: Optional[Dict] = None,
        mutation_interpreter: Optional[MutationInterpreter] = None
    ):
        self.mutation = mutation
        self.pdb_id = pdb_id
        self.chain = chain
        self.position = position
        self.wt_aa = wt_aa
        self.mut_aa = mut_aa
        self.foldx = foldx_results
        self.alphamissense = alphamissense_score
        self.structure = structure_context or {}
        self.interpreter = mutation_interpreter
        
        # Get amino acid property comparison
        self.aa_comparison = compare_amino_acids(wt_aa, mut_aa)
    
    def generate_complete_report(self) -> Dict:
        """
        Generate comprehensive mutation analysis report
        
        Returns:
            Complete report with all sections
        """
        return {
            'metadata': self._generate_metadata(),
            'location_context': self._generate_location_context(),
            'predictions': self._generate_predictions(),
            'structural_impact': self._generate_structural_impact(),
            'energetic_analysis': self._generate_energetic_analysis(),
            'biochemical_analysis': self._generate_biochemical_analysis(),
            'interpretations': self._generate_interpretations(),
            'summary': self._generate_summary(),
            'clinical_assessment': self._generate_clinical_assessment()
        }
     
    def _generate_metadata(self) -> Dict:
        """Basic mutation information"""
        return {
            'mutation': self.mutation,
            'pdb_id': self.pdb_id,
            'chain': self.chain,
            'position': self.position,
            'wild_type': {
                'code': self.wt_aa,
                'name': AMINO_ACIDS[self.wt_aa]['name'],
                'three_letter': AMINO_ACIDS[self.wt_aa]['three_letter']
            },
            'mutant': {
                'code': self.mut_aa,
                'name': AMINO_ACIDS[self.mut_aa]['name'],
                'three_letter': AMINO_ACIDS[self.mut_aa]['three_letter']
            }
        }
    
    def _generate_location_context(self) -> Dict:
        """Structural location information"""
        if not self.structure:
            return {'available': False}
        
        location = self.structure.get('location', {})
        ss = self.structure.get('secondary_structure', {})
        
        return {
            'available': True,
            'location': location.get('location', 'unknown'),
            'description': location.get('description', ''),
            'sasa': location.get('sasa', None),
            'secondary_structure': {
                'type': ss.get('type', 'unknown'),
                'code': ss.get('code', ''),
                'description': ss.get('description', '')
            },
            'flexibility': self.structure.get('bfactor', {})
        }
    
    def _generate_predictions(self) -> Dict:
        """Prediction scores and interpretations"""
        predictions = {}
        
        # FoldX stability prediction
        if self.foldx:
            ddg = self.foldx.get('ddg', 0.0)
            sd = self.foldx.get('sd', 0.0)
            
            predictions['stability'] = {
                'source': 'FoldX',
                'ddg': ddg,
                'sd': sd,
                'ddg_formatted': f"{ddg:+.2f} ± {sd:.2f} kcal/mol",
                'confidence': 'high' if sd < 0.5 else 'medium' if sd < 1.0 else 'low',
                'num_runs': self.foldx.get('num_runs', 1),
                'interpretation': self._interpret_foldx_ddg(ddg),
                'category': self._categorize_ddg(ddg)
            }
        
        # AlphaMissense pathogenicity prediction
        if self.alphamissense is not None:
            predictions['pathogenicity'] = {
                'source': 'AlphaMissense',
                'score': self.alphamissense,
                'score_formatted': f"{self.alphamissense:.3f}",
                'interpretation': self._interpret_alphamissense(self.alphamissense),
                'category': self._categorize_alphamissense(self.alphamissense)
            }
        
        return predictions
     
    def _generate_structural_impact(self) -> Dict:
        """Property-based structural analysis"""
        diff = self.aa_comparison['differences']
        
        return {
            'size_change': {
                'value': diff['size']['value'],
                'description': diff['size']['description'],
                'percentage': (diff['size']['value'] / 5.0) * 100  # Out of max 5
            },
            'charge_change': {
                'changed': diff['charge']['changed'],
                'type': diff['charge']['type'],
                'from': self.aa_comparison['reference']['charge'],
                'to': self.aa_comparison['alternative']['charge']
            },
            'hydrophobicity_change': {
                'value': diff['hydrophobicity']['value'],
                'description': diff['hydrophobicity']['description'],
                'from': self.aa_comparison['reference']['hydrophobicity'],
                'to': self.aa_comparison['alternative']['hydrophobicity']
            },
            'polarity_change': {
                'changed': diff['polarity']['changed'],
                'from': self.aa_comparison['reference']['polarity'],
                'to': self.aa_comparison['alternative']['polarity']
            },
            'special_properties': {
                'aromatic_change': diff['aromatic_change'],
                'sulphur_change': diff['sulphur_change'],
                'hydroxyl_change': diff['hydroxyl_change']
            }
        }
     
    def _generate_energetic_analysis(self) -> Dict:
        """FoldX energy breakdown analysis"""
        if not self.foldx or 'energy_breakdown' not in self.foldx:
            return {'available': False}
        
        breakdown = self.foldx['energy_breakdown']
        
        # Identify main contributors
        contributors = []
        for term, value in breakdown.items():
            if abs(value) > 0.05:  # Significant contribution threshold
                contributors.append({
                    'term': term,
                    'value': value,
                    'value_formatted': f"{value:+.2f} kcal/mol",
                    'type': 'destabilizing' if value > 0 else 'stabilizing',
                    'interpretation': self._interpret_energy_term(term, value)
                })
        
        # Sort by absolute magnitude
        contributors.sort(key=lambda x: abs(x['value']), reverse=True)
        
        # Separate stabilizing and destabilizing
        stabilizing = [c for c in contributors if c['value'] < 0]
        destabilizing = [c for c in contributors if c['value'] > 0]
        
        return {
            'available': True,
            'total_ddg': self.foldx.get('ddg', 0.0),
            'all_terms': breakdown,
            'main_contributors': contributors[:5],  # Top 5
            'stabilizing_factors': stabilizing,
            'destabilizing_factors': destabilizing,
            'net_structural': sum(v for k, v in breakdown.items() if 'entropy' not in k.lower()),
            'net_entropy': sum(v for k, v in breakdown.items() if 'entropy' in k.lower())
        }
       
    def _generate_biochemical_analysis(self) -> Dict:
        """HOPE-style biochemical interpretations"""
        interpretations = []
        
        # Use MutationInterpreter if available
        if self.interpreter:
            full_analysis = self.interpreter.analyze_mutation(
                wt_residue=self.wt_aa,
                mut_residue=self.mut_aa,
                chain=self.chain,
                position=self.position,
                alphamissense_score=self.alphamissense,
                foldx_ddg=self.foldx.get('ddg') if self.foldx else None
            )
            
            # Convert to our format
            for interp in full_analysis.get('interpretations', []):
                interpretations.append(Interpretation(
                    category=interp['category'],
                    severity=Severity(interp['severity']),
                    icon=interp['icon'],
                    text=interp['text']
                ))
        else:
            # Fallback: basic analysis
            interpretations = self._generate_basic_interpretations()
        
        return {
            'interpretations': [
                {
                    'category': i.category,
                    'severity': i.severity.value,
                    'icon': i.icon,
                    'text': i.text
                }
                for i in interpretations
            ],
            'property_changes': self.aa_comparison['differences']
        }
      
    def _generate_interpretations(self) -> List[Dict]:
        """Combined interpretation list"""
        return self._generate_biochemical_analysis()['interpretations']
       
    def _generate_summary(self) -> str:
        """Executive summary of mutation impact"""
        parts = []
        
        # Overall assessment
        severity_count = self._count_severity_levels()
        if severity_count['high'] >= 2:
            parts.append("This mutation is predicted to be highly damaging to protein structure and function.")
        elif severity_count['high'] == 1 or severity_count['medium'] >= 2:
            parts.append("This mutation is likely to impact protein structure or function.")
        else:
            parts.append("This mutation may have mild to moderate effects on protein structure.")
        
        # FoldX assessment
        if self.foldx:
            ddg = self.foldx.get('ddg', 0.0)
            if ddg > 2.0:
                parts.append(f"FoldX predicts significant destabilization (ΔΔG: +{ddg:.2f} kcal/mol).")
            elif ddg > 1.0:
                parts.append(f"FoldX predicts moderate destabilization (ΔΔG: +{ddg:.2f} kcal/mol).")
            elif ddg < -1.0:
                parts.append(f"FoldX predicts stabilization (ΔΔG: {ddg:.2f} kcal/mol).")
            else:
                parts.append(f"FoldX predicts minimal stability change (ΔΔG: {ddg:.2f} kcal/mol).")
        
        # AlphaMissense assessment
        if self.alphamissense is not None:
            if self.alphamissense > 0.7:
                parts.append(f"AlphaMissense predicts this as pathogenic (score: {self.alphamissense:.3f}).")
            elif self.alphamissense < 0.3:
                parts.append(f"AlphaMissense predicts this as benign (score: {self.alphamissense:.3f}).")
            else:
                parts.append(f"AlphaMissense gives an ambiguous prediction (score: {self.alphamissense:.3f}).")
        
        # Location context
        if self.structure:
            location = self.structure.get('location', {}).get('location')
            if location == 'core':
                parts.append("The mutation occurs in the protein core, which increases the likelihood of structural disruption.")
            elif location == 'surface':
                parts.append("The mutation is on the protein surface, which may reduce structural impact.")
        
        return " ".join(parts)
      
    def _generate_clinical_assessment(self) -> Dict:
        """Overall clinical significance assessment"""
        # Combine all evidence
        evidence_score = 0
        confidence_factors = []
        
        # FoldX contribution
        if self.foldx:
            ddg = self.foldx.get('ddg', 0.0)
            if ddg > 2.0:
                evidence_score += 3
                confidence_factors.append("High destabilization (FoldX)")
            elif ddg > 1.0:
                evidence_score += 2
                confidence_factors.append("Moderate destabilization (FoldX)")
            elif ddg < -1.0:
                evidence_score -= 2
                confidence_factors.append("Stabilizing effect (FoldX)")
        
        # AlphaMissense contribution
        if self.alphamissense is not None:
            if self.alphamissense > 0.7:
                evidence_score += 3
                confidence_factors.append("Pathogenic prediction (AlphaMissense)")
            elif self.alphamissense < 0.3:
                evidence_score -= 2
                confidence_factors.append("Benign prediction (AlphaMissense)")
        
        # Severity contribution
        severity_count = self._count_severity_levels()
        evidence_score += severity_count['high'] * 2
        evidence_score += severity_count['medium'] * 1
        
        # Determine overall classification
        if evidence_score >= 5:
            classification = "Likely Pathogenic"
            confidence = "High"
        elif evidence_score >= 3:
            classification = "Possibly Pathogenic"
            confidence = "Medium"
        elif evidence_score <= -3:
            classification = "Likely Benign"
            confidence = "High"
        elif evidence_score <= -1:
            classification = "Possibly Benign"
            confidence = "Medium"
        else:
            classification = "Uncertain Significance"
            confidence = "Low"
        
        return {
            'classification': classification,
            'confidence': confidence,
            'evidence_score': evidence_score,
            'supporting_evidence': confidence_factors,
            'recommendation': self._generate_recommendation(classification, confidence)
        }
     
    # ==================== Helper Methods ====================
    
    def _interpret_foldx_ddg(self, ddg: float) -> str:
        """Human-readable FoldX interpretation"""
        if ddg > 2.0:
            return "Highly destabilizing. This mutation significantly disrupts protein stability."
        elif ddg > 1.0:
            return "Moderately destabilizing. This mutation may affect protein stability."
        elif ddg > -1.0:
            return "Neutral effect. Minimal impact on protein stability."
        else:
            return "Stabilizing. This mutation enhances protein stability."
    
    def _categorize_ddg(self, ddg: float) -> str:
        """Categorize ΔΔG value"""
        if ddg > 2.0:
            return "highly_destabilizing"
        elif ddg > 1.0:
            return "moderately_destabilizing"
        elif ddg > -1.0:
            return "neutral"
        else:
            return "stabilizing"
    
    def _interpret_alphamissense(self, score: float) -> str:
        """Human-readable AlphaMissense interpretation"""
        if score > 0.7:
            return f"High probability of pathogenicity. This mutation is likely disease-causing."
        elif score > 0.5:
            return f"Moderate probability of pathogenicity. This mutation may be disease-related."
        elif score > 0.3:
            return f"Ambiguous prediction. Additional evidence needed for classification."
        else:
            return f"Low probability of pathogenicity. This mutation is likely benign."
    
    def _categorize_alphamissense(self, score: float) -> str:
        """Categorize AlphaMissense score"""
        if score > 0.7:
            return "pathogenic"
        elif score > 0.5:
            return "likely_pathogenic"
        elif score > 0.3:
            return "uncertain"
        else:
            return "benign"
    
    def _interpret_energy_term(self, term: str, value: float) -> str:
        """Interpret individual energy term change"""
        interpretations = {
            'backbone_hbond': {
                'pos': "Loss of backbone hydrogen bonds reduces stability",
                'neg': "Gain of backbone hydrogen bonds enhances stability"
            },
            'sidechain_hbond': {
                'pos': "Loss of sidechain hydrogen bonds is destabilizing",
                'neg': "Formation of new sidechain hydrogen bonds"
            },
            'van_der_waals': {
                'pos': "Loss of van der Waals contacts weakens structure",
                'neg': "Improved van der Waals packing strengthens structure"
            },
            'electrostatics': {
                'pos': "Unfavorable electrostatic interactions introduced",
                'neg': "Favorable electrostatic interactions formed"
            },
            'solvation_polar': {
                'pos': "Increased polar solvation penalty (burial of polar groups)",
                'neg': "Improved polar solvation (exposure of polar groups)"
            },
            'solvation_hydrophobic': {
                'pos': "Poor hydrophobic burial increases energy",
                'neg': "Better hydrophobic burial lowers energy"
            },
            'clash': {
                'pos': "Steric clashes introduced between atoms",
                'neg': "Reduced steric clashes"
            },
            'entropy_sidechain': {
                'pos': "Loss of sidechain conformational freedom",
                'neg': "Increased sidechain conformational freedom"
            },
            'entropy_mainchain': {
                'pos': "Loss of backbone flexibility",
                'neg': "Increased backbone flexibility"
            }
        }
        
        term_interp = interpretations.get(term, {})
        if value > 0:
            return term_interp.get('pos', f"{term} increased by {value:.2f}")
        else:
            return term_interp.get('neg', f"{term} decreased by {value:.2f}")
     
    def _generate_basic_interpretations(self) -> List[Interpretation]:
        """Generate basic interpretations without MutationInterpreter"""
        interpretations = []
        diff = self.aa_comparison['differences']
        
        # Charge analysis
        if diff['charge']['changed']:
            if diff['charge']['type'] == 'charge_reversal':
                interpretations.append(Interpretation(
                    category='charge',
                    severity=Severity.HIGH,
                    icon='⚠️',
                    text=f"Charge reversal ({self.aa_comparison['reference']['charge']} → "
                         f"{self.aa_comparison['alternative']['charge']}) likely disrupts "
                         f"electrostatic interactions."
                ))
        
        # Size analysis
        size_diff = diff['size']['value']
        if abs(size_diff) >= 2:
            if size_diff > 0:
                interpretations.append(Interpretation(
                    category='size',
                    severity=Severity.MEDIUM,
                    icon='📏',
                    text=f"Larger residue may cause steric clashes."
                ))
            else:
                interpretations.append(Interpretation(
                    category='size',
                    severity=Severity.MEDIUM,
                    icon='🕳️',
                    text=f"Smaller residue may create destabilizing cavity."
                ))
        
        return interpretations
    
    def _count_severity_levels(self) -> Dict[str, int]:
        """Count interpretations by severity"""
        interps = self._generate_biochemical_analysis()['interpretations']
        return {
            'high': sum(1 for i in interps if i['severity'] == 'high'),
            'medium': sum(1 for i in interps if i['severity'] == 'medium'),
            'low': sum(1 for i in interps if i['severity'] == 'low')
        }
    
    def _generate_recommendation(self, classification: str, confidence: str) -> str:
        """Generate clinical recommendation"""
        recommendations = {
            ('Likely Pathogenic', 'High'): 
                "Strong evidence for pathogenicity. Recommend experimental validation and clinical correlation.",
            ('Likely Pathogenic', 'Medium'):
                "Substantial evidence for pathogenicity. Consider additional functional studies.",
            ('Possibly Pathogenic', 'Medium'):
                "Moderate evidence for pathogenicity. Further investigation recommended.",
            ('Uncertain Significance', 'Low'):
                "Insufficient evidence for classification. Requires additional functional and population data.",
            ('Likely Benign', 'High'):
                "Strong evidence for benign effect. Mutation likely tolerated.",
            ('Likely Benign', 'Medium'):
                "Substantial evidence for benign effect. Likely tolerated but monitor for edge cases.",
        }
        
        return recommendations.get(
            (classification, confidence),
            "Classification uncertain. Additional evidence required."
        )
     
    # ==================== Export Methods ====================
    
    def to_json(self) -> Dict:
        """Export complete report as JSON"""
        return self.generate_complete_report()
    
    def to_markdown(self) -> str:
        """Export report as formatted markdown"""
        report = self.generate_complete_report()
        
        md = f"# Mutation Analysis: {self.mutation}\n\n"
        md += f"**{self.aa_comparison['reference']['name']} → "
        md += f"{self.aa_comparison['alternative']['name']}** "
        md += f"at Chain {self.chain}, Position {self.position}\n\n"
        
        md += "---\n\n"
        
        # Location
        if report['location_context']['available']:
            loc = report['location_context']
            md += f"## 📍 Location Context\n\n"
            md += f"- **Location:** {loc['description']}\n"
            if loc['sasa']:
                md += f"- **SASA:** {loc['sasa']:.1f} Ų\n"
            md += f"- **Secondary Structure:** {loc['secondary_structure']['description']}\n\n"
        
        # Predictions
        md += "## 🔬 Predictions\n\n"
        
        if 'pathogenicity' in report['predictions']:
            path = report['predictions']['pathogenicity']
            md += f"### AlphaMissense\n"
            md += f"- **Score:** {path['score_formatted']}\n"
            md += f"- **Classification:** {path['category']}\n"
            md += f"- {path['interpretation']}\n\n"
        
        if 'stability' in report['predictions']:
            stab = report['predictions']['stability']
            md += f"### FoldX Stability\n"
            md += f"- **ΔΔG:** {stab['ddg_formatted']}\n"
            md += f"- **Confidence:** {stab['confidence']}\n"
            md += f"- {stab['interpretation']}\n\n"
        
        # Summary
        md += f"## 💡 Summary\n\n"
        md += f"{report['summary']}\n\n"
        
        # Clinical assessment
        md += f"## ✅ Clinical Assessment\n\n"
        clin = report['clinical_assessment']
        md += f"**Classification:** {clin['classification']}\n\n"
        md += f"**Confidence:** {clin['confidence']}\n\n"
        md += f"**Recommendation:** {clin['recommendation']}\n\n"
        
        return md

# ==================== Usage Example ====================
"""
if __name__ == "__main__":
    # Example usage with your FoldX results
    
    foldx_results = {
        'ddg': 0.230,
        'sd': 0.0,
        'wt_energy': 133.243,
        'mut_energy': 133.473,
        'num_runs': 5,
        'energy_breakdown': {
            'backbone_hbond': 0.0,
            'sidechain_hbond': 0.0,
            'van_der_waals': -0.024,
            'electrostatics': 0.0,
            'solvation_polar': -0.076,
            'solvation_hydrophobic': -0.186,
            'clash': 0.005,
            'entropy_sidechain': 0.198,
            'entropy_mainchain': 0.312
        }
    }
    
    structure_context = {
        'location': {
            'location': 'surface',
            'sasa': 45.2,
            'description': 'exposed on protein surface'
        },
        'secondary_structure': {
            'type': 'helix',
            'code': 'H',
            'description': 'alpha helix'
        }
    }
    
    # Create report
    report = DataLensReport(
        mutation="SA729A",
        pdb_id="7lmk",
        chain="A",
        position=729,
        wt_aa="S",
        mut_aa="A",
        foldx_results=foldx_results,
        alphamissense_score=0.12,
        structure_context=structure_context
    )
    
    # Generate complete report
    full_report = report.generate_complete_report()
    
    # Export as JSON
    import json
    print(json.dumps(full_report, indent=2))
    
    # Export as Markdown
    print("\n" + "="*60 + "\n")
    print(report.to_markdown())

    """