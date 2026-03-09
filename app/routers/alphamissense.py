"""
app/routers/alphamissense.py
AlphaMissense prediction endpoints
Uses alphamissense_service to query alphamisense_predictions table
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Optional

from app.services.alphamissense_service import (
    check_uniprot_has_predictions,
    get_all_variants_for_pdb,
    get_prediction_for_mutation,  # This matches the service function name
    get_all_predictions_at_position,
    get_chains_with_am_data,
    get_prediction_count
)

router = APIRouter(prefix="/api/alphamissense", tags=["alphamissense"])


class CoverageCheck(BaseModel):
    """AlphaMissense coverage status"""
    uniprot_id: str
    has_predictions: bool
    prediction_count: int


class MutationPrediction(BaseModel):
    """AlphaMissense prediction for a mutation"""
    uniprot_id: str
    protein_variant: str
    position: int
    reference_aa: str
    alternate_aa: str
    am_pathogenicity: float
    am_class: str


class PositionPredictions(BaseModel):
    """All predictions at a position"""
    uniprot_id: str
    position: int
    reference_aa: str
    predictions: List[Dict]
    prediction_count: int


class AnalyzableChains(BaseModel):
    """Chains that have AlphaMissense data"""
    pdb_id: str
    chains: List[Dict]
    chain_count: int


@router.get("/coverage/{uniprot_id}", response_model=CoverageCheck)
async def check_coverage(uniprot_id: str):
    """
    Check if AlphaMissense predictions exist for a UniProt ID
    
    Example: GET /api/alphamissense/coverage/P62805
    Returns: {"uniprot_id": "P62805", "has_predictions": true, "prediction_count": 1957}
    
    Example: GET /api/alphamissense/coverage/Q24K09
    Returns: {"uniprot_id": "Q24K09", "has_predictions": false, "prediction_count": 0}
    
    Use this to check BEFORE querying for predictions!
    """
    has_predictions = await check_uniprot_has_predictions(uniprot_id)
    count = await get_prediction_count(uniprot_id) if has_predictions else 0
    
    return CoverageCheck(
        uniprot_id=uniprot_id.upper(),
        has_predictions=has_predictions,
        prediction_count=count
    )


@router.get("/predict/{uniprot_id}/{position}/{alternate_aa}")
async def predict_mutation(uniprot_id: str, position: int, alternate_aa: str):
    """
    Get AlphaMissense prediction for a specific mutation
    
    Example: GET /api/alphamissense/predict/P62805/28/H
    Returns pathogenicity score and classification
    
    Returns 404 if no prediction found (either no coverage or mutation not in DB)
    """
    # First check if this UniProt has ANY predictions
    has_coverage = await check_uniprot_has_predictions(uniprot_id)
    
    if not has_coverage:
        raise HTTPException(
            status_code=404,
            detail=f"No AlphaMissense predictions for UniProt {uniprot_id}. "
                   f"This protein is not in the AlphaMissense dataset."
        )
    
    # Get specific prediction
    prediction = await get_prediction_for_mutation(uniprot_id, position, alternate_aa)
    
    if not prediction:
        raise HTTPException(
            status_code=404,
            detail=f"No prediction found for {uniprot_id} position {position} variant {alternate_aa}"
        )
    
    return MutationPrediction(**prediction)


@router.get("/position/{uniprot_id}/{position}", response_model=PositionPredictions)
async def get_all_mutations_at_position(uniprot_id: str, position: int):
    """
    Get all possible mutations at a position (all 19 amino acid substitutions)
    Sorted by pathogenicity score (most pathogenic first)
    
    Example: GET /api/alphamissense/position/P62805/28
    Returns all mutations: K28A, K28C, K28D, ..., K28Y with scores
    """
    # Check coverage first
    has_coverage = await check_uniprot_has_predictions(uniprot_id)
    
    if not has_coverage:
        raise HTTPException(
            status_code=404,
            detail=f"No AlphaMissense predictions for UniProt {uniprot_id}"
        )
    
    predictions = await get_all_predictions_at_position(uniprot_id, position)
    
    if not predictions:
        raise HTTPException(
            status_code=404,
            detail=f"No predictions found for position {position}"
        )
    
    reference_aa = predictions[0]['reference_aa'] if predictions else '?'
    
    return PositionPredictions(
        uniprot_id=uniprot_id.upper(),
        position=position,
        reference_aa=reference_aa,
        predictions=predictions,
        prediction_count=len(predictions)
    )


@router.get("/pdb/{pdb_id}/analyzable-chains", response_model=AnalyzableChains)
async def get_chains_with_predictions(pdb_id: str):

    """
    Get ONLY chains that have AlphaMissense prediction data
    
    Example: GET /api/alphamissense/pdb/7lmk/analyzable-chains
    
    For 7lmk:
    - Returns chains A,B (P62805) ✓ Has predictions
    - Excludes chains C,D (Q24K09) ✗ No predictions
    
    Use this to automatically select the best chain for analysis!
    """
    chains = await get_chains_with_am_data(pdb_id)
    
    if not chains:
        raise HTTPException(
            status_code=404,
            detail=f"No chains with AlphaMissense data for PDB {pdb_id}. "
                   f"This structure may only contain proteins not in AlphaMissense dataset."
        )
    
    return AnalyzableChains(
        pdb_id=pdb_id.upper(),
        chains=chains,
        chain_count=len(chains)
    )


@router.get("/variants/{pdb_id}")    # ← no /api prefix
async def get_all_variants_for_pdb_route(pdb_id: str):
    """
    Get all AlphaMissense variants for a PDB ID via JOIN query.
    Used to populate the variant selector panel on PDB load.
    
    Example: GET /api/alphamissense/variants/7lmk
    """
    rows = await get_all_variants_for_pdb(pdb_id)
    
    if not rows:
        raise HTTPException(
            status_code=404,
            detail=f"No variant data found for PDB {pdb_id}"
        )
    
    return {"variants": rows}