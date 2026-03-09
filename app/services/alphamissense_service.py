"""
app/services/alphamissense_service.py
AlphaMissense prediction queries with validation
Uses table: alphamisense_predictions
"""

from typing import Optional, Dict, List
# db_pool = None
from app.database import init_db_pool


async def check_uniprot_has_predictions(uniprot_id: str) -> bool:
    """Check if UniProt ID has predictions (P62805=True, Q24K09=False)"""
    pool = await init_db_pool()
    query = "SELECT EXISTS(SELECT 1 FROM alphamisense_predictions WHERE uniprot_id = $1 LIMIT 1)"
    async with pool.acquire() as conn:
        return await conn.fetchval(query, uniprot_id)


async def get_prediction_count(uniprot_id: str) -> int:
    """Get total number of predictions for a UniProt ID"""
    pool = await init_db_pool()
    query = "SELECT COUNT(*) FROM alphamisense_predictions WHERE uniprot_id = $1"
    async with pool.acquire() as conn:
        return await conn.fetchval(query, uniprot_id)


async def get_prediction_for_mutation(uniprot_id: str, position: int, alternate_aa: str) -> Optional[Dict]:
    """Get AlphaMissense prediction for specific mutation"""
    pool = await init_db_pool()
    
    query = """
        SELECT 
            uniprot_id, protein_variant, am_pathogenicity, am_class,
            position, reference_aa, alternate_aa
        FROM alphamisense_predictions
        WHERE uniprot_id = $1 AND position = $2 AND alternate_aa = $3
        LIMIT 1
    """
    
    async with pool.acquire() as conn:
        result = await conn.fetchrow(query, uniprot_id, position, alternate_aa.upper())
        
        if result:
            return {
                'uniprot_id': result['uniprot_id'],
                'protein_variant': result['protein_variant'],
                'am_pathogenicity': float(result['am_pathogenicity']),
                'am_class': result['am_class'],
                'position': result['position'],
                'reference_aa': result['reference_aa'],
                'alternate_aa': result['alternate_aa']
            }
        return None


async def get_all_predictions_at_position(uniprot_id: str, position: int) -> List[Dict]:
    """Get all 19 possible substitutions at a position"""
    pool = await init_db_pool()
    
    query = """
        SELECT protein_variant, am_pathogenicity, am_class, reference_aa, alternate_aa
        FROM alphamisense_predictions
        WHERE uniprot_id = $1 AND position = $2
        ORDER BY am_pathogenicity DESC
    """
    
    async with pool.acquire() as conn:
        rows = await conn.fetch(query, uniprot_id, position)
        
        return [{
            'protein_variant': row['protein_variant'],
            'am_pathogenicity': float(row['am_pathogenicity']),
            'am_class': row['am_class'],
            'reference_aa': row['reference_aa'],
            'alternate_aa': row['alternate_aa']
        } for row in rows]


async def get_chains_with_am_data(pdb_id: str) -> List[Dict]:
    """
    CRITICAL: Get ONLY chains that have AlphaMissense data
    INNER JOIN filters out Q24K09 automatically
    """
    pool = await init_db_pool()
    
    query = """
        SELECT DISTINCT
            m.chain_id, m.uniprot_id, m.pdb_start, m.pdb_end,
            m.uniprot_start, m.uniprot_end,
            COUNT(a.protein_variant) as prediction_count
        FROM uniprot_pdb_mapping m
        INNER JOIN alphamisense_predictions a ON m.uniprot_id = a.uniprot_id
        WHERE m.pdb_id = $1
        GROUP BY m.chain_id, m.uniprot_id, m.pdb_start, m.pdb_end, m.uniprot_start, m.uniprot_end
        ORDER BY m.chain_id
    """
    
    async with pool.acquire() as conn:
        rows = await conn.fetch(query, pdb_id.lower())
        return [dict(row) for row in rows]

async def get_all_variants_for_pdb(pdb_id: str):
    pool = await init_db_pool()   # ← same as every other function

    query = """
    SELECT DISTINCT
        ap.protein_variant  AS variant_id,
        ap.position,
        ap.reference_aa     AS ref_aa,
        ap.alternate_aa     AS alt_aa,
        ap.am_pathogenicity AS pathogenicity_score,
        ap.am_class         AS classification
    FROM uniprot_pdb_mapping pum
    JOIN alphamisense_predictions ap ON ap.uniprot_id = pum.uniprot_id
    WHERE pum.pdb_id = lower($1)
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