# app/routers/structure.py
from fastapi import APIRouter, HTTPException
from app.models import (
    FetchPDBRequest, 
    FetchPDBResponse, 
    StructureMetadata
)
from app.services.structure_service import StructureService
from app.utils.logger import logger

# Create router - THIS LINE IS CRITICAL
router = APIRouter(prefix="/structure", tags=["structure_router"])

# Initialize service
structure_service = StructureService()

@router.post("/fetch", response_model=FetchPDBResponse)
async def fetch_pdb(request: FetchPDBRequest):
    """
    Fetch PDB structure from RCSB, repair if needed
    """
    try:
        pdb_content, metadata = await structure_service.fetch_and_prepare(
            request.pdb_id,
            force_repair=request.force_repair
        )
        
        struct_metadata = StructureMetadata(**metadata)
        
        return FetchPDBResponse(
            success=True,
            pdb_content=pdb_content,
            metadata=struct_metadata,
            message=f"Structure {request.pdb_id} prepared successfully"
        )
    
    except ValueError as e:
        logger.error(f"ValueError fetching {request.pdb_id}: {e}")
        raise HTTPException(status_code=404, detail=str(e))
    
    except Exception as e:
        logger.error(f"Unexpected error fetching {request.pdb_id}: {e}")
        raise HTTPException(status_code=500, detail="Internal server error")

@router.get("/cache/info")
async def cache_info():
    """Get cache statistics"""
    from app.config import settings
    
    cache_files = list(settings.CACHE_DIR.glob("*.pdb"))
    total_size = sum(f.stat().st_size for f in cache_files)
    
    return {
        "cached_structures": len(cache_files),
        "total_size_mb": total_size / (1024 * 1024),
        "cache_dir": str(settings.CACHE_DIR)
    }