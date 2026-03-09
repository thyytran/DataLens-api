# app/services/structure_service.py
import urllib.request
from typing import Tuple, Optional
from app.config import settings
from app.services.foldx_service import FoldXService
from app.utils.logger import logger
import time

class StructureService:
    def __init__(self):
        self.foldx = FoldXService()
        self.cache_dir = settings.CACHE_DIR
    
    async def fetch_and_prepare(
        self, 
        pdb_id: str, 
        force_repair: bool = False
    ) -> Tuple[str, dict]:
        """
        Complete pipeline: Download → Repair → Return
        """
        pdb_id = pdb_id.upper()
        logger.info(f"Processing PDB {pdb_id}")
        
        # Check cache
        if not force_repair and settings.CACHE_ENABLED:
            cached = self._get_cached(pdb_id)
            if cached:
                logger.info(f"Using cached structure for {pdb_id}")
                return cached
        
        # Download original
        logger.info(f"Downloading {pdb_id} from RCSB")
        original_pdb = await self._download_pdb(pdb_id)
        
        if not original_pdb:
            raise ValueError(f"Failed to download PDB {pdb_id}")
        
        # Always repair (like MutationExplorer does for consistency)
        logger.info(f"Running FoldX RepairPDB for {pdb_id}")
        
        start_time = time.time()
        
        repaired_pdb, energy = await self.foldx.repair_pdb(original_pdb, pdb_id)
        
        repair_time = time.time() - start_time
        
        # Count atoms/residues
        atom_count = original_pdb.count('\nATOM')
        residue_count = len(self._extract_residues(original_pdb))
        
        # Cache the result
        if settings.CACHE_ENABLED:
            self._cache_structure(pdb_id, repaired_pdb, energy)
        
        metadata = {
            "pdb_id": pdb_id,
            "source": "repaired",
            "energy": energy,
            "repair_time_seconds": repair_time,
            "atom_count": atom_count,
            "residue_count": residue_count
        }
        
        logger.info(f"Successfully prepared {pdb_id} in {repair_time:.2f}s")
        
        return repaired_pdb, metadata
    
    async def _download_pdb(self, pdb_id: str) -> Optional[str]:
        """Download PDB from RCSB"""
        url = f"{settings.RCSB_PDB_URL}/{pdb_id}.pdb"
        
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                content = response.read().decode('utf-8')
                logger.info(f"Downloaded {pdb_id}: {len(content)} bytes")
                return content
        except Exception as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            return None
    
    def _get_cached(self, pdb_id: str) -> Optional[Tuple[str, dict]]:
        """Retrieve from cache if exists and not expired"""
        cached_file = self.cache_dir / f"{pdb_id}_Repair.pdb"
        metadata_file = self.cache_dir / f"{pdb_id}_metadata.json"
        
        if not cached_file.exists():
            return None
        
        # Check age
        import time
        age_days = (time.time() - cached_file.stat().st_mtime) / (24 * 3600)
        if age_days > settings.CACHE_MAX_AGE_DAYS:
            logger.info(f"Cache expired for {pdb_id} ({age_days:.1f} days old)")
            return None
        
        # Read cached structure
        with open(cached_file, 'r') as f:
            pdb_content = f.read()
        
        # Read metadata if exists
        metadata = {"pdb_id": pdb_id, "source": "cache"}
        if metadata_file.exists():
            import json
            with open(metadata_file, 'r') as f:
                metadata.update(json.load(f))
        
        return pdb_content, metadata
    
    def _cache_structure(self, pdb_id: str, pdb_content: str, energy: dict):
        """Save structure to cache"""
        cached_file = self.cache_dir / f"{pdb_id}_Repair.pdb"
        metadata_file = self.cache_dir / f"{pdb_id}_metadata.json"
        
        # Save PDB
        with open(cached_file, 'w') as f:
            f.write(pdb_content)
        
        # Save metadata
        import json
        metadata = {
            "pdb_id": pdb_id,
            "energy": energy,
            "cached_at": time.time()
        }
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        logger.info(f"Cached {pdb_id}")
    
    def _extract_residues(self, pdb_content: str) -> set:
        """Extract unique residues from PDB"""
        residues = set()
        for line in pdb_content.split('\n'):
            if line.startswith('ATOM'):
                res_id = line[17:27]  # residue name + number + chain
                residues.add(res_id)
        return residues
    
    