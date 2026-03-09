# backend/routers/mutation_analysis.py

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import Optional

from app.services import foldx_service, structure_service

# Only import what we actually use
#from services.mutation_interpreter import MutationInterpreter
#from services.datalens_reports import DataLensReport
# ← Removed FoldXResult - it's not used!

router = APIRouter(prefix="/api", tags=["mutation_analysis"])
foldx_service = foldx_service.FoldXService()
structure_service = structure_service.StructureService()


class MutationRequest(BaseModel):
    pdb_id: str
    chain: str
    position: int
    wt_aa: str
    mut_aa: str
    alphamissense_score: Optional[float] = None


@router.post("/analyze_mutation")
async def analyze_mutation(request: MutationRequest):
    """
    Full FoldX pipeline:
    1. Download + RepairPDB
    2. BuildModel with selected mutation
    3. Return ΔΔG + summary
    """
    pdb_id = request.pdb_id.upper()

    # ── Step 1: Get repaired PDB (from cache or run RepairPDB) ──────────────
    try:
        repaired_pdb, metadata = await structure_service.fetch_and_prepare(pdb_id)
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"RepairPDB failed: {str(e)}")

    # ── Step 2: BuildModel with selected mutation ────────────────────────────
    # FoldX mutation format: "WtAAChainPositionMutAA;"
    # e.g. L A 60 G → "LA60G;"
    foldx_mutation = f"{request.wt_aa}{request.chain}{request.position}{request.mut_aa};"


    try:
        mutant_pdb, ddg = await foldx_service.build_model(
            pdb_content=repaired_pdb,
            pdb_id=pdb_id,
            mutations=[foldx_mutation]
        )
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"BuildModel failed: {str(e)}")

    # ── Step 3: Interpret result ─────────────────────────────────────────────
    interpretation = _interpret_ddg(ddg)
    am_interpretation = _interpret_am(request.alphamissense_score)

    return {
        "mutation":             foldx_mutation.rstrip(";"),
        "pdb_id":               pdb_id,
        "chain":                request.chain,
        "position":             request.position,
        "wt_aa":                request.wt_aa,
        "mut_aa":               request.mut_aa,
        "mutant_pdb": mutant_pdb,   # ← add this, it's already returned by build_mode
        "predictions": {
            "foldx": {
                "ddg":            ddg,
                "interpretation": interpretation
            },
            "alphamissense": {
                "score":          request.alphamissense_score,
                "interpretation": am_interpretation
            }
        },
        "summary": (
            f"{request.wt_aa}{request.position}{request.mut_aa}: "
            f"ΔΔG = {ddg:+.2f} kcal/mol ({interpretation}), "
            f"AlphaMissense = {request.alphamissense_score:.3f} ({am_interpretation})"
            if request.alphamissense_score else
            f"{request.wt_aa}{request.position}{request.mut_aa}: "
            f"ΔΔG = {ddg:+.2f} kcal/mol ({interpretation})"
        )
    }


def _interpret_ddg(ddg: float) -> str:
    """Interpret FoldX ΔΔG value"""
    if ddg > 2.0:
        return "highly destabilizing"
    elif ddg > 0.5:
        return "destabilizing"
    elif ddg > -0.5:
        return "neutral"
    elif ddg > -2.0:
        return "stabilizing"
    else:
        return "highly stabilizing"


def _interpret_am(score: Optional[float]) -> str:
    """Interpret AlphaMissense score"""
    if score is None:
        return "no data"
    if score >= 0.564:
        return "likely pathogenic"
    elif score >= 0.340:
        return "ambiguous"
    else:
        return "likely benign"