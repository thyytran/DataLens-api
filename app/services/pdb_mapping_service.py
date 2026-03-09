"""
app/services/pdb_mapping_service.py
Query PDB to UniProt mappings from PostgreSQL database
Uses existing uniprot_pdb_mapping table
"""

from typing import Dict, Optional, List
from app.database import init_db_pool


async def get_pdb_chain_mapping(pdb_id: str, chain_id: str) -> Optional[Dict]:
    """
    Get UniProt ID and mapping details for a specific PDB chain
    
    Args:
        pdb_id: PDB identifier (e.g., '7lmk')
        chain_id: Chain identifier (e.g., 'A')
    
    Returns:
        Dictionary with mapping details or None if not found
    """
    pool = await init_db_pool()
    
    query = """
        SELECT 
            pdb_id, 
            chain_id, 
            uniprot_id,
            pdb_start,
            pdb_end,
            observed_start,
            observed_end,
            uniprot_start,
            uniprot_end
        FROM uniprot_pdb_mapping 
        WHERE pdb_id = $1 AND chain_id = $2
        LIMIT 1
    """
    
    async with pool.acquire() as conn:
        result = await conn.fetchrow(query, pdb_id.lower(), chain_id.upper())
        
        if result:
            return {
                'pdb_id': result['pdb_id'],
                'chain_id': result['chain_id'],
                'uniprot_id': result['uniprot_id'],
                'pdb_start': result['pdb_start'],
                'pdb_end': result['pdb_end'],
                'observed_start': result['observed_start'],
                'observed_end': result['observed_end'],
                'uniprot_start': result['uniprot_start'],
                'uniprot_end': result['uniprot_end']
            }
        
        return None

async def get_all_chains_for_pdb(pdb_id: str) -> List[Dict]:
    """
    Get all chain → UniProt mappings for a PDB entry with details
    
    Args:
        pdb_id: PDB identifier (e.g., '7lmk')
    
    Returns:
        List of dictionaries with mapping details for each chain
    """
    pool = await init_db_pool()
    
    query = """
        SELECT 
            chain_id, 
            uniprot_id,
            pdb_start,
            pdb_end,
            observed_start,
            observed_end,
            uniprot_start,
            uniprot_end
        FROM uniprot_pdb_mapping 
        WHERE pdb_id = $1
        ORDER BY chain_id
    """
    
    async with pool.acquire() as conn:
        rows = await conn.fetch(query, pdb_id.lower())
        
        mappings = []
        for row in rows:
            mappings.append({
                'chain_id': row['chain_id'],
                'uniprot_id': row['uniprot_id'],
                'pdb_start': row['pdb_start'],
                'pdb_end': row['pdb_end'],
                'observed_start': row['observed_start'],
                'observed_end': row['observed_end'],
                'uniprot_start': row['uniprot_start'],
                'uniprot_end': row['uniprot_end']
            })
        
        return mappings

async def map_pdb_to_uniprot_residue(pdb_id: str, chain_id: str, pdb_residue: int) -> Optional[int]:
    """
    Convert PDB residue number to UniProt residue number
    
    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        pdb_residue: Residue number in PDB coordinates
    
    Returns:
        UniProt residue number or None if out of range
    """
    mapping = await get_pdb_chain_mapping(pdb_id, chain_id)
    
    if not mapping:
        return None
    
    # Check if PDB residue is in the mapped range
    if mapping['pdb_start'] <= pdb_residue <= mapping['pdb_end']:
        # Calculate offset
        offset = pdb_residue - mapping['pdb_start']
        uniprot_residue = mapping['uniprot_start'] + offset
        return uniprot_residue
    
    return None

async def map_uniprot_to_pdb_residue(pdb_id: str, chain_id: str, uniprot_residue: int) -> Optional[int]:
    """
    Convert UniProt residue number to PDB residue number
    
    Args:
        pdb_id: PDB identifier
        chain_id: Chain identifier
        uniprot_residue: Residue number in UniProt coordinates
    
    Returns:
        PDB residue number or None if out of range
    """
    mapping = await get_pdb_chain_mapping(pdb_id, chain_id)
    
    if not mapping:
        return None
    
    # Check if UniProt residue is in the mapped range
    if mapping['uniprot_start'] <= uniprot_residue <= mapping['uniprot_end']:
        # Calculate offset
        offset = uniprot_residue - mapping['uniprot_start']
        pdb_residue = mapping['pdb_start'] + offset
        return pdb_residue
    
    return None

async def get_all_chains_with_residue_ranges(pdb_id: str) -> List[Dict]:
    """
    Get all chains and their residue ranges from uniprot_pdb_mapping.
    No AM filter — includes all UniProt IDs mapped to this PDB.
    Used for Section 2 of variant selector (FoldX manual mutation).
    """
    pool = await init_db_pool()

    query = """
        SELECT
            pdb_id,
            chain_id,
            uniprot_id,
            uniprot_start AS pdb_start,
            uniprot_end AS pdb_end,
            observed_start,
            observed_end,
            pdb_start AS uniprot_start,
            pdb_end AS uniprot_end
        FROM uniprot_pdb_mapping
        WHERE pdb_id = lower($1)
        ORDER BY chain_id, pdb_start
    """

    async with pool.acquire() as conn:
        rows = await conn.fetch(query, pdb_id)
        return [dict(r) for r in rows]

async def close_db_pool():
    """Close database connection pool"""
    global db_pool
    if db_pool:
        await db_pool.close()
        db_pool = None