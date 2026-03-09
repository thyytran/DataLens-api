# amino_acid_properties.py

AMINO_ACIDS = {
    # Single-letter code: {properties}
    'A': {
        'name': 'Alanine',
        'three_letter': 'ALA',
        'molecular_weight': 89.09,
        'charge': 'neutral',
        'polarity': 'nonpolar',
        'hydrophobicity': 1.8,  # Kyte-Doolittle scale
        'size': 'tiny',
        'size_value': 1,  # Relative scale 1-5
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'C': {
        'name': 'Cysteine',
        'three_letter': 'CYS',
        'molecular_weight': 121.15,
        'charge': 'neutral',
        'polarity': 'polar',
        'hydrophobicity': 2.5,
        'size': 'small',
        'size_value': 2,
        'aromatic': False,
        'sulphur': True,
        'hydroxyl': False
    },
    'D': {
        'name': 'Aspartic acid',
        'three_letter': 'ASP',
        'molecular_weight': 133.10,
        'charge': 'negative',
        'polarity': 'polar',
        'hydrophobicity': -3.5,
        'size': 'small',
        'size_value': 2,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'E': {
        'name': 'Glutamic acid',
        'three_letter': 'GLU',
        'molecular_weight': 147.13,
        'charge': 'negative',
        'polarity': 'polar',
        'hydrophobicity': -3.5,
        'size': 'medium',
        'size_value': 3,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'F': {
        'name': 'Phenylalanine',
        'three_letter': 'PHE',
        'molecular_weight': 165.19,
        'charge': 'neutral',
        'polarity': 'nonpolar',
        'hydrophobicity': 2.8,
        'size': 'large',
        'size_value': 4,
        'aromatic': True,
        'sulphur': False,
        'hydroxyl': False
    },
    'G': {
        'name': 'Glycine',
        'three_letter': 'GLY',
        'molecular_weight': 75.07,
        'charge': 'neutral',
        'polarity': 'nonpolar',
        'hydrophobicity': -0.4,
        'size': 'tiny',
        'size_value': 1,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'H': {
        'name': 'Histidine',
        'three_letter': 'HIS',
        'molecular_weight': 155.16,
        'charge': 'positive',  # Can be positive at physiological pH
        'polarity': 'polar',
        'hydrophobicity': -3.2,
        'size': 'medium',
        'size_value': 3,
        'aromatic': True,
        'sulphur': False,
        'hydroxyl': False
    },
    'I': {
        'name': 'Isoleucine',
        'three_letter': 'ILE',
        'molecular_weight': 131.17,
        'charge': 'neutral',
        'polarity': 'nonpolar',
        'hydrophobicity': 4.5,
        'size': 'large',
        'size_value': 4,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'K': {
        'name': 'Lysine',
        'three_letter': 'LYS',
        'molecular_weight': 146.19,
        'charge': 'positive',
        'polarity': 'polar',
        'hydrophobicity': -3.9,
        'size': 'large',
        'size_value': 4,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'L': {
        'name': 'Leucine',
        'three_letter': 'LEU',
        'molecular_weight': 131.17,
        'charge': 'neutral',
        'polarity': 'nonpolar',
        'hydrophobicity': 3.8,
        'size': 'large',
        'size_value': 4,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'M': {
        'name': 'Methionine',
        'three_letter': 'MET',
        'molecular_weight': 149.21,
        'charge': 'neutral',
        'polarity': 'nonpolar',
        'hydrophobicity': 1.9,
        'size': 'large',
        'size_value': 4,
        'aromatic': False,
        'sulphur': True,
        'hydroxyl': False
    },
    'N': {
        'name': 'Asparagine',
        'three_letter': 'ASN',
        'molecular_weight': 132.12,
        'charge': 'neutral',
        'polarity': 'polar',
        'hydrophobicity': -3.5,
        'size': 'small',
        'size_value': 2,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'P': {
        'name': 'Proline',
        'three_letter': 'PRO',
        'molecular_weight': 115.13,
        'charge': 'neutral',
        'polarity': 'nonpolar',
        'hydrophobicity': -1.6,
        'size': 'small',
        'size_value': 2,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'Q': {
        'name': 'Glutamine',
        'three_letter': 'GLN',
        'molecular_weight': 146.15,
        'charge': 'neutral',
        'polarity': 'polar',
        'hydrophobicity': -3.5,
        'size': 'medium',
        'size_value': 3,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'R': {
        'name': 'Arginine',
        'three_letter': 'ARG',
        'molecular_weight': 174.20,
        'charge': 'positive',
        'polarity': 'polar',
        'hydrophobicity': -4.5,
        'size': 'very_large',
        'size_value': 5,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'S': {
        'name': 'Serine',
        'three_letter': 'SER',
        'molecular_weight': 105.09,
        'charge': 'neutral',
        'polarity': 'polar',
        'hydrophobicity': -0.8,
        'size': 'tiny',
        'size_value': 1,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': True
    },
    'T': {
        'name': 'Threonine',
        'three_letter': 'THR',
        'molecular_weight': 119.12,
        'charge': 'neutral',
        'polarity': 'polar',
        'hydrophobicity': -0.7,
        'size': 'small',
        'size_value': 2,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': True
    },
    'V': {
        'name': 'Valine',
        'three_letter': 'VAL',
        'molecular_weight': 117.15,
        'charge': 'neutral',
        'polarity': 'nonpolar',
        'hydrophobicity': 4.2,
        'size': 'medium',
        'size_value': 3,
        'aromatic': False,
        'sulphur': False,
        'hydroxyl': False
    },
    'W': {
        'name': 'Tryptophan',
        'three_letter': 'TRP',
        'molecular_weight': 204.23,
        'charge': 'neutral',
        'polarity': 'nonpolar',
        'hydrophobicity': -0.9,
        'size': 'very_large',
        'size_value': 5,
        'aromatic': True,
        'sulphur': False,
        'hydroxyl': False
    },
    'Y': {
        'name': 'Tyrosine',
        'three_letter': 'TYR',
        'molecular_weight': 181.19,
        'charge': 'neutral',
        'polarity': 'polar',
        'hydrophobicity': -1.3,
        'size': 'large',
        'size_value': 4,
        'aromatic': True,
        'sulphur': False,
        'hydroxyl': True
    }
}


# Mutation notation parser
def parse_mutation(mutation_string):
    """
    Parse mutation notation like 'L60G' or 'LEU60GLY'
    
    Returns:
        dict: {
            'reference_aa': 'L',
            'position': 60,
            'alternative_aa': 'G',
            'reference_name': 'Leucine',
            'alternative_name': 'Glycine'
        }
    """
    import re
    
    # Handle single-letter notation (L60G)
    match = re.match(r'([A-Z])(\d+)([A-Z])', mutation_string)
    if match:
        ref_aa = match.group(1)
        position = int(match.group(2))
        alt_aa = match.group(3)
        
        return {
            'reference_aa': ref_aa,
            'position': position,
            'alternative_aa': alt_aa,
            'reference_name': AMINO_ACIDS[ref_aa]['name'],
            'alternative_name': AMINO_ACIDS[alt_aa]['name'],
            'reference_3letter': AMINO_ACIDS[ref_aa]['three_letter'],
            'alternative_3letter': AMINO_ACIDS[alt_aa]['three_letter']
        }
    
    # Handle three-letter notation (LEU60GLY)
    match = re.match(r'([A-Z]{3})(\d+)([A-Z]{3})', mutation_string)
    if match:
        ref_3letter = match.group(1)
        position = int(match.group(2))
        alt_3letter = match.group(3)
        
        # Find single-letter codes
        ref_aa = next(k for k, v in AMINO_ACIDS.items() if v['three_letter'] == ref_3letter)
        alt_aa = next(k for k, v in AMINO_ACIDS.items() if v['three_letter'] == alt_3letter)
        
        return {
            'reference_aa': ref_aa,
            'position': position,
            'alternative_aa': alt_aa,
            'reference_name': AMINO_ACIDS[ref_aa]['name'],
            'alternative_name': AMINO_ACIDS[alt_aa]['name'],
            'reference_3letter': ref_3letter,
            'alternative_3letter': alt_3letter
        }
    
    raise ValueError(f"Invalid mutation notation: {mutation_string}")


# Property comparison function
def compare_amino_acids(ref_aa, alt_aa):
    """
    Compare properties between reference and alternative amino acids
    
    Returns:
        dict: Comprehensive comparison of all properties
    """
    ref_props = AMINO_ACIDS[ref_aa]
    alt_props = AMINO_ACIDS[alt_aa]
    
    # Calculate differences
    weight_diff = alt_props['molecular_weight'] - ref_props['molecular_weight']
    size_diff = alt_props['size_value'] - ref_props['size_value']
    hydrophobicity_diff = alt_props['hydrophobicity'] - ref_props['hydrophobicity']
    
    # Determine charge change
    charge_change = None
    if ref_props['charge'] != alt_props['charge']:
        if ref_props['charge'] == 'neutral':
            charge_change = 'gain_charge'
        elif alt_props['charge'] == 'neutral':
            charge_change = 'loss_charge'
        else:
            charge_change = 'charge_reversal'  # pos->neg or neg->pos
    
    # Determine polarity change
    polarity_change = ref_props['polarity'] != alt_props['polarity']
    
    return {
        'reference': ref_props,
        'alternative': alt_props,
        'differences': {
            'molecular_weight': {
                'value': weight_diff,
                'percentage': (weight_diff / ref_props['molecular_weight']) * 100
            },
            'size': {
                'value': size_diff,
                'description': 'larger' if size_diff > 0 else 'smaller' if size_diff < 0 else 'same'
            },
            'hydrophobicity': {
                'value': hydrophobicity_diff,
                'description': 'more_hydrophobic' if hydrophobicity_diff > 0 
                              else 'more_hydrophilic' if hydrophobicity_diff < 0 
                              else 'similar'
            },
            'charge': {
                'changed': charge_change is not None,
                'type': charge_change
            },
            'polarity': {
                'changed': polarity_change,
                'from': ref_props['polarity'],
                'to': alt_props['polarity']
            },
            'aromatic_change': ref_props['aromatic'] != alt_props['aromatic'],
            'sulphur_change': ref_props['sulphur'] != alt_props['sulphur'],
            'hydroxyl_change': ref_props['hydroxyl'] != alt_props['hydroxyl']
        }
    }

