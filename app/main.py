# app/main.py
from typing import Dict
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.config import settings
from app.utils.logger import logger
from app.models import HealthCheckResponse
from app.database import init_db_pool

from app.routers.alphamissense import router as alphamissense_router
from app.routers.pdb_mapping import router as pdb_mapping_router
from app.routers.uniprot import router as uniprot_router
from app.routers.mutation_analysis import router as mutation_analysis_router
from app.routers.structure import router as structure_router

db_pool = None
# In-memory cache for protein data
protein_cache: Dict[str, dict] = {}

app = FastAPI(
    title=settings.APP_NAME,
    description="Backend service for DataLens protein mutation analysis",
    version="1.0.0"
)

app.include_router(alphamissense_router)
app.include_router(pdb_mapping_router)
app.include_router(uniprot_router)
app.include_router(mutation_analysis_router)
app.include_router(structure_router)

# CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Import router AFTER app is created
from app.routers.structure import router as structure_router

# Include router
app.include_router(structure_router)

@app.on_event("startup")
async def startup_event():
    await init_db_pool()
    logger.info(f"Starting {settings.APP_NAME}")
    logger.info(f"FoldX path: {settings.FOLDX_PATH}")
    logger.info(f"Cache directory: {settings.CACHE_DIR}")

@app.on_event("shutdown")
async def shutdown_event():
    logger.info("Shutting down")

@app.get("/", tags=["health"])
async def root():
    return {
        "app": settings.APP_NAME,
        "status": "running",
        "version": "1.0.0"
    }

@app.get("/health", response_model=HealthCheckResponse, tags=["health"])
async def health_check():
    """Health check endpoint"""
    from pathlib import Path
    
    # Check FoldX
    foldx_available = Path(settings.FOLDX_PATH).exists()
    
    # Check cache
    cache_size = 0
    try:
        cache_size = sum(
            f.stat().st_size for f in settings.CACHE_DIR.glob("*") if f.is_file()
        ) / (1024 * 1024)  # MB
    except:
        pass
    
    return HealthCheckResponse(
        status="healthy" if foldx_available else "degraded",
        foldx_available=foldx_available,
        database_connected=True,
        cache_size_mb=cache_size
    )
