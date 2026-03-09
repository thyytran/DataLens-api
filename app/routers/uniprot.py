"""
app/routers/uniprot.py
UniProt API endpoints
Uses uniprot_service to fetch and parse UniProt REST API data
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Optional

from app.services.uniprot_service import (
    fetch_uniprot_data,
    parse_secondary_structure,
    parse_variants,
    parse_ptms,
    parse_dna_binding,
    get_protein_name,
    get_function,
    get_sequence,
    get_secondary_structure_for_position,
    is_dna_binding_position,
    clear_cache,
    get_cache_status
)

router = APIRouter(prefix="/api/uniprot", tags=["uniprot"])


class ResidueInfo(BaseModel):
    """Residue information from UniProt"""
    position: int
    amino_acid: str
    secondary_structure: Dict[str, str]
    dna_binding: bool
    variants: List[Dict]
    ptms: List[Dict]


class ProteinInfo(BaseModel):
    """Full protein information from UniProt"""
    uniprot_id: str
    protein_name: Optional[str]
    sequence: Optional[str]
    sequence_length: int
    function: Optional[str]
    secondary_structure: List[Dict]
    variants: List[Dict]
    ptms: List[Dict]
    dna_binding: List[Dict]


@router.get("/protein/{uniprot_id}", response_model=ProteinInfo)
async def get_protein_info(uniprot_id: str):
    """
    Fetch full protein information from UniProt
    Caches result for subsequent queries
    
    Example: GET /api/uniprot/protein/P62805
    Returns complete protein data including function, variants, PTMs, etc.
    """
    try:
        data = fetch_uniprot_data(uniprot_id)
        sequence = get_sequence(data)
        
        return ProteinInfo(
            uniprot_id=uniprot_id.upper(),
            protein_name=get_protein_name(data),
            sequence=sequence,
            sequence_length=len(sequence) if sequence else 0,
            function=get_function(data),
            secondary_structure=parse_secondary_structure(data),
            variants=parse_variants(data),
            ptms=parse_ptms(data),
            dna_binding=parse_dna_binding(data)
        )
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to fetch UniProt data: {str(e)}"
        )


@router.get("/protein/{uniprot_id}/residue/{position}", response_model=ResidueInfo)
async def get_residue_info(uniprot_id: str, position: int):
    """
    Get comprehensive information for a specific residue
    
    Example: GET /api/uniprot/protein/P62805/residue/28
    Returns secondary structure, variants, PTMs for residue 28
    
    Note: Position should be in UniProt coordinates!
    """
    try:
        data = fetch_uniprot_data(uniprot_id)
        
        # Get sequence and validate position
        sequence = get_sequence(data)
        if not sequence or position < 1 or position > len(sequence):
            raise HTTPException(
                status_code=400,
                detail=f"Invalid position {position} "
                       f"(sequence length: {len(sequence) if sequence else 0})"
            )
        
        # Parse all features
        ss_list = parse_secondary_structure(data)
        all_variants = parse_variants(data)
        all_ptms = parse_ptms(data)
        binding_list = parse_dna_binding(data)
        
        # Get residue-specific info
        amino_acid = sequence[position - 1]
        secondary_structure = get_secondary_structure_for_position(ss_list, position)
        dna_binding = is_dna_binding_position(binding_list, position)
        
        # Filter variants and PTMs for this position
        variants = [v for v in all_variants if v['position'] == position]
        ptms = [p for p in all_ptms if p['position'] == position]
        
        return ResidueInfo(
            position=position,
            amino_acid=amino_acid,
            secondary_structure=secondary_structure,
            dna_binding=dna_binding,
            variants=variants,
            ptms=ptms[:5]  # Limit to first 5 PTMs
        )
    
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Failed to process residue info: {str(e)}"
        )


@router.get("/cache/status")
async def get_cache_status_endpoint():
    """
    Get UniProt cache statistics
    
    Shows which proteins are currently cached in memory
    """
    return get_cache_status()


@router.post("/cache/clear")
async def clear_cache_endpoint():
    """
    Clear UniProt cache
    
    Forces fresh fetch from UniProt API on next request
    """
    count = clear_cache()
    return {"message": f"Cleared {count} cached proteins"}