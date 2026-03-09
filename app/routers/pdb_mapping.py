"""
app/routers/pdb_mapping.py
PDB to UniProt mapping endpoints
Uses pdb_mapping_service to query uniprot_pdb_mapping table
"""

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Dict, Optional

from app.services.pdb_mapping_service import (
    get_all_chains_with_residue_ranges,
    get_pdb_chain_mapping,
    get_all_chains_for_pdb,
    map_pdb_to_uniprot_residue,
    map_uniprot_to_pdb_residue
)

router = APIRouter(prefix="/api", tags=["pdb-mapping"])


class ChainMapping(BaseModel):
    """Chain mapping details"""
    pdb_id: str
    chain_id: str
    uniprot_id: str
    pdb_start: int
    pdb_end: int
    observed_start: int
    observed_end: int
    uniprot_start: int
    uniprot_end: int


class PDBChainList(BaseModel):
    """List of chains for a PDB"""
    pdb_id: str
    chains: List[ChainMapping]


class ResidueMapping(BaseModel):
    """PDB to UniProt residue conversion"""
    pdb_id: str
    chain_id: str
    pdb_residue: int
    uniprot_residue: int
    uniprot_id: str


@router.get("/pdb/{pdb_id}/chains", response_model=PDBChainList)
async def get_pdb_chains(pdb_id: str):
    """
    Get all chains for a PDB entry
    
    Example: GET /api/pdb/7lmk/chains
    Returns all chains with their UniProt mappings and residue ranges
    """
    chains = await get_all_chains_for_pdb(pdb_id)
    
    if not chains:
        raise HTTPException(
            status_code=404,
            detail=f"PDB {pdb_id} not found in mapping database"
        )
    
    # Convert list of dicts to ChainMapping objects
    chain_mappings = []
    for chain in chains:
        chain_mappings.append(ChainMapping(
            pdb_id=pdb_id.lower(),
            **chain
        ))
    
    return PDBChainList(
        pdb_id=pdb_id.lower(),
        chains=chain_mappings
    )


@router.get("/pdb/{pdb_id}/chain/{chain_id}", response_model=ChainMapping)
async def get_chain_detail(pdb_id: str, chain_id: str):
    """
    Get detailed mapping for a specific chain
    
    Example: GET /api/pdb/7lmk/chain/A
    Returns UniProt ID and residue number mappings for chain A
    """
    mapping = await get_pdb_chain_mapping(pdb_id, chain_id)
    
    if not mapping:
        raise HTTPException(
            status_code=404,
            detail=f"No mapping found for PDB {pdb_id} chain {chain_id}"
        )
    
    return ChainMapping(
        pdb_id=pdb_id.lower(),
        **mapping
    )


@router.get("/pdb/{pdb_id}/chain/{chain_id}/residue/{pdb_residue}", response_model=ResidueMapping)
async def convert_residue_number(pdb_id: str, chain_id: str, pdb_residue: int):
    """
    Convert PDB residue number to UniProt residue number
    
    Example: GET /api/pdb/7lmk/chain/A/residue/60
    Returns: {"pdb_residue": 60, "uniprot_residue": 28, "uniprot_id": "P62805"}
    
    Critical for coordinate mapping!
    """
    # Get mapping details
    mapping = await get_pdb_chain_mapping(pdb_id, chain_id)
    
    if not mapping:
        raise HTTPException(
            status_code=404,
            detail=f"No mapping for PDB {pdb_id} chain {chain_id}"
        )
    
    # Convert coordinate
    uniprot_residue = await map_pdb_to_uniprot_residue(pdb_id, chain_id, pdb_residue)
    
    if uniprot_residue is None:
        raise HTTPException(
            status_code=400,
            detail=f"PDB residue {pdb_residue} is out of mapped range "
                   f"({mapping['pdb_start']}-{mapping['pdb_end']})"
        )
    
    return ResidueMapping(
        pdb_id=pdb_id.lower(),
        chain_id=chain_id,
        pdb_residue=pdb_residue,
        uniprot_residue=uniprot_residue,
        uniprot_id=mapping['uniprot_id']
    )

@router.get("/pdb/{pdb_id}/all-residues")
async def get_all_residues(pdb_id: str):
    """
    Get all chains and residue ranges for a PDB from uniprot_pdb_mapping.
    Includes ALL UniProt IDs — not filtered by AlphaMissense coverage.
    Used to populate Section 2 of the variant selector for FoldX mutations.

    Example: GET /api/pdb/7lmk/all-residues
    Returns chains A,B (P62805) AND chains C,D (Q24K09)
    """
    chains = await get_all_chains_with_residue_ranges(pdb_id)

    if not chains:
        raise HTTPException(
            status_code=404,
            detail=f"No chain mapping found for PDB {pdb_id}"
        )

    return {
        "pdb_id":      pdb_id.lower(),
        "chain_count": len(chains),
        "chains":      chains
    }