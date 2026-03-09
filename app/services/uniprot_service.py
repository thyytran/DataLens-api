"""
app/services/uniprot_service.py
UniProt data fetching and parsing service for DataLens
"""

from typing import Optional, List, Dict
import requests
from ..config import settings

# Constants
UNIPROT_BASE_URL = settings.UNIPROT_BASE_URL

# In-memory cache for protein data
_protein_cache: Dict[str, dict] = {}


def fetch_uniprot_data(uniprot_id: str) -> dict:
    """Fetch data from UniProt REST API with caching"""
    uniprot_id = uniprot_id.upper()
    
    # Check cache first
    if uniprot_id in _protein_cache:
        return _protein_cache[uniprot_id]
    
    # Fetch from API
    url = f"{UNIPROT_BASE_URL}/{uniprot_id}.json"
    
    try:
        response = requests.get(url, timeout=10)
        
        if response.status_code == 200:
            data = response.json()
            _protein_cache[uniprot_id] = data
            return data
        else:
            raise Exception(f"UniProt API returned status {response.status_code}")
    
    except requests.exceptions.Timeout:
        raise Exception("UniProt API timeout")
    except requests.exceptions.RequestException as e:
        raise Exception(f"Network error: {str(e)}")

def parse_secondary_structure(data: dict) -> List[Dict]:
    """Extract secondary structure annotations"""
    if 'features' not in data:
        return []
    
    ss_list = []
    for feature in data['features']:
        if feature['type'] in ['Helix', 'Beta strand', 'Turn']:
            ss_list.append({
                'type': feature['type'],
                'start': feature['location']['start']['value'],
                'end': feature['location']['end']['value']
            })
    
    return ss_list

def parse_variants(data: dict) -> List[Dict]:
    """Extract known disease/somatic variants"""
    if 'features' not in data:
        return []
    
    variant_list = []
    for feature in data['features']:
        if feature['type'] == 'Natural variant':
            try:
                description = feature.get('description', 'Unknown')
                variant = {
                    'position': feature['location']['start']['value'],
                    'from_aa': feature['alternativeSequence']['originalSequence'],
                    'to_aa': feature['alternativeSequence']['alternativeSequences'][0],
                    'description': description,
                    'disease': _extract_disease_name(description)
                }
                variant_list.append(variant)
            except (KeyError, IndexError):
                continue
    
    return variant_list

def parse_ptms(data: dict) -> List[Dict]:
    """Extract post-translational modifications"""
    if 'features' not in data:
        return []
    
    ptm_list = []
    for feature in data['features']:
        if feature['type'] == 'Modified residue':
            try:
                ptm = {
                    'position': feature['location']['start']['value'],
                    'modification': feature['description']
                }
                ptm_list.append(ptm)
            except KeyError:
                continue
    
    return ptm_list

def parse_dna_binding(data: dict) -> List[Dict]:
    """Extract DNA binding regions"""
    if 'features' not in data:
        return []
    
    binding_list = []
    for feature in data['features']:
        if feature['type'] == 'DNA binding':
            try:
                binding = {
                    'start': feature['location']['start']['value'],
                    'end': feature['location']['end']['value']
                }
                binding_list.append(binding)
            except KeyError:
                continue
    
    return binding_list

def get_protein_name(data: dict) -> Optional[str]:
    """Extract recommended protein name"""
    try:
        return data['proteinDescription']['recommendedName']['fullName']['value']
    except KeyError:
        return None

def get_function(data: dict) -> Optional[str]:
    """Extract protein function description"""
    if 'comments' not in data:
        return None
    
    for comment in data['comments']:
        if comment['commentType'] == 'FUNCTION':
            try:
                return comment['texts'][0]['value']
            except (KeyError, IndexError):
                continue
    
    return None

def get_sequence(data: dict) -> Optional[str]:
    """Extract protein sequence"""
    try:
        return data['sequence']['value']
    except KeyError:
        return None

def get_secondary_structure_for_position(ss_list: List[Dict], position: int) -> Dict[str, str]:
    """Get secondary structure type for a specific position"""
    for ss in ss_list:
        if ss['start'] <= position <= ss['end']:
            return {
                'type': ss['type'],
                'range': f"{ss['start']}-{ss['end']}"
            }
    
    return {'type': 'Loop/Coil', 'range': ''}

def is_dna_binding_position(binding_list: List[Dict], position: int) -> bool:
    """Check if position is in DNA binding region"""
    for binding in binding_list:
        if binding['start'] <= position <= binding['end']:
            return True
    return False

def _extract_disease_name(description: str) -> str:
    """Extract disease name from variant description"""
    if 'in ' in description:
        parts = description.split('in ')
        if len(parts) > 1:
            disease = parts[1].split(';')[0].split('.')[0]
            return disease.strip()
    return 'Unknown'

def clear_cache():
    """Clear the protein data cache"""
    count = len(_protein_cache)
    _protein_cache.clear()
    return count

def get_cache_status():
    """Get cache statistics"""
    return {
        "cached_proteins": len(_protein_cache),
        "uniprot_ids": list(_protein_cache.keys())
    }