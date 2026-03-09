# app/models.py
from pydantic import BaseModel, Field
from typing import Optional, List, Dict
from enum import Enum

class StructureSource(str, Enum):
    ORIGINAL = "original"
    REPAIRED = "repaired"
    CACHE = "cache"

class FetchPDBRequest(BaseModel):
    pdb_id: str = Field(..., min_length=4, max_length=4, description="PDB ID (4 characters)")
    force_repair: bool = Field(False, description="Force repair even if cached")

class StructureMetadata(BaseModel):
    pdb_id: str
    source: StructureSource
    issues_found: Optional[List[str]] = None
    energy: Optional[Dict[str, float]] = None
    repair_time_seconds: Optional[float] = None
    atom_count: int
    residue_count: int

class FetchPDBResponse(BaseModel):
    success: bool
    pdb_content: str
    metadata: StructureMetadata
    message: Optional[str] = None

class MutationRequest(BaseModel):
    pdb_id: str
    pdb_content: Optional[str] = None  # If already loaded
    chain: str
    position: int
    wildtype: str = Field(..., min_length=1, max_length=3)
    mutant: str = Field(..., min_length=1, max_length=3)

class MutationResult(BaseModel):
    success: bool
    foldx_ddg: Optional[float] = None
    alphamisense_score: Optional[float] = None
    prediction: Optional[str] = None  # "stabilizing", "destabilizing", "neutral"
    message: Optional[str] = None

class HealthCheckResponse(BaseModel):
    status: str
    foldx_available: bool
    database_connected: bool
    cache_size_mb: float

# Response models
class ChainMapping(BaseModel):
    chain_id: str
    uniprot_id: str
    pdb_start: str
    pdb_end: str
    uniprot_start: int
    uniprot_end: int
    coverage_length: int

class ChainSummary(BaseModel):
    chain_id: str
    uniprot_id: str
    uniprot_start: int
    uniprot_end: int
    positions_with_data: int
    total_mutations: int
    avg_pathogenicity: float
    pathogenic_count: int
    benign_count: int
    ambiguous_count: int

class PDBStructureData(BaseModel):
    pdb_id: str
    chains: List[ChainMapping]
    chain_summaries: List[ChainSummary]
    total_mutations: int

class ResidueInfo(BaseModel):
    """Residue information response model"""
    position: int
    amino_acid: str
    secondary_structure: Dict[str, str]
    dna_binding: bool
    variants: List[Dict]
    ptms: List[Dict]


class ProteinInfo(BaseModel):
    """Full protein information response model"""
    uniprot_id: str
    protein_name: Optional[str]
    sequence: Optional[str]
    sequence_length: int
    function: Optional[str]
    secondary_structure: List[Dict]
    variants: List[Dict]
    ptms: List[Dict]
    dna_binding: List[Dict]