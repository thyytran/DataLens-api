# app/services/foldx_service.py
import subprocess
import tempfile
import asyncio
from pathlib import Path
from typing import List, Tuple, Dict
from app.config import settings
from app.utils.logger import logger

class FoldXService:
    def __init__(self):
        self.foldx_path = settings.FOLDX_PATH
        #self.rotabase_path = settings.ROTABASE_PATH
        
        # Validate FoldX installation
        if not Path(self.foldx_path).exists():
            raise FileNotFoundError(f"FoldX not found at {self.foldx_path}")
        
        #if not Path(self.rotabase_path).exists():
        #    raise FileNotFoundError(f"rotabase.txt not found at {self.rotabase_path}")
    
    async def repair_pdb(self, pdb_content: str, pdb_id: str) -> Tuple[str, Dict]:
        """
        Run FoldX RepairPDB asynchronously
        
        Returns:
            (repaired_pdb_content, energy_dict)
        """
        # Use a thread pool for blocking subprocess call
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(
            None,
            self._repair_pdb_sync,
            pdb_content,
            pdb_id
        )
    
    def _repair_pdb_sync(self, pdb_content: str, pdb_id: str) -> Tuple[str, Dict]:
        """Synchronous FoldX RepairPDB execution"""
        
        with tempfile.TemporaryDirectory(dir=settings.TEMP_DIR) as tmpdir:
            tmpdir = Path(tmpdir)
            
            # Write input PDB
            input_pdb = tmpdir / f"{pdb_id}.pdb"
            with open(input_pdb, 'w') as f:
                f.write(pdb_content)
            
            # Copy rotabase
            #shutil.copy(self.rotabase_path, tmpdir / "rotabase.txt")
            
            # Build command
            cmd = [
                str(self.foldx_path),
                "--command=RepairPDB",
                f"--pdb={pdb_id}.pdb",
                "--ionStrength=0.05",
                "--temperature=298",
                "--pH=7.0"
            ]
            
            logger.debug(f"Running: {' '.join(cmd)}")
            
            try:
                # Run FoldX
                result = subprocess.run(
                    cmd,
                    cwd=tmpdir,
                    capture_output=True,
                    text=True,
                    timeout=settings.FOLDX_TIMEOUT_SECONDS
                )
                
                # Check output
                repaired_pdb = tmpdir / f"{pdb_id}_Repair.pdb"
                
                if not repaired_pdb.exists():
                    logger.error(f"FoldX failed for {pdb_id}")
                    logger.error(f"STDOUT: {result.stdout}")
                    logger.error(f"STDERR: {result.stderr}")
                    raise RuntimeError(f"FoldX RepairPDB failed for {pdb_id}")
                
                # Read repaired structure
                with open(repaired_pdb, 'r') as f:
                    repaired_content = f.read()
                
                # Parse energy
                energy = self._parse_foldx_output(tmpdir / f"Repair_{pdb_id}.fxout")
                
                logger.info(f"RepairPDB successful for {pdb_id}: {energy}")
                
                return repaired_content, energy
            
            except subprocess.TimeoutExpired:
                logger.error(f"FoldX timeout for {pdb_id}")
                raise RuntimeError(f"FoldX timeout after {settings.FOLDX_TIMEOUT_SECONDS}s")
            
            except Exception as e:
                logger.error(f"FoldX error for {pdb_id}: {e}")
                raise
    
    def _parse_foldx_output(self, fxout_file: Path) -> Dict[str, float]:
        """Parse FoldX energy output"""
        if not fxout_file.exists():
            return {}
        
        energy = {}
        with open(fxout_file, 'r') as f:
            for line in f:
                if 'Total energy' in line:
                    parts = line.strip().split()
                    try:
                        energy['total'] = float(parts[-1])
                    except (ValueError, IndexError):
                        pass
                
                # Parse individual energy terms if needed
                # backbone_hbond, sidechain_hbond, van_der_waals, etc.
        
        return energy
    
    async def build_model(
        self, 
        pdb_content: str, 
        pdb_id: str,
        mutations: List[str]
    ) -> Tuple[str, float]:
        """
        Run FoldX BuildModel for mutations
        
        Args:
            pdb_content: Repaired PDB content
            pdb_id: PDB identifier
            mutations: List of mutations in format "RA620H;" (chain+residue+new)
        
        Returns:
            (mutant_pdb_content, ddg)
        """
        loop = asyncio.get_event_loop()
        return await loop.run_in_executor(
            None,
            self._build_model_sync,
            pdb_content,
            pdb_id,
            mutations
        )
    
    def _build_model_sync(
        self,
        pdb_content: str,
        pdb_id: str,
        mutations: List[str]
    ) -> Tuple[str, float]:
        with tempfile.TemporaryDirectory(dir=settings.TEMP_DIR) as tmpdir:
            tmpdir = Path(tmpdir)

            # Write repaired PDB
            input_pdb = tmpdir / f"{pdb_id}_Repair.pdb"
            with open(input_pdb, 'w') as f:
                f.write(pdb_content)

            # Write mutation list inside tmpdir
            mutation_file = tmpdir / "individual_list.txt"
            with open(mutation_file, 'w') as f:
                for mut in mutations:
                    f.write(mut + '\n')

            cmd = [
                str(self.foldx_path),
                "--command=BuildModel",
                f"--pdb={pdb_id}_Repair.pdb",
                "--mutant-file=individual_list.txt",   # relative to cwd
                "--numberOfRuns=1"
            ]

            logger.debug(f"Running BuildModel: {' '.join(cmd)}")

            try:
                result = subprocess.run(
                    cmd,
                    cwd=tmpdir,          # FoldX runs here, finds both files
                    capture_output=True,
                    text=True,
                    timeout=settings.FOLDX_TIMEOUT_SECONDS
                )

                logger.debug(f"FoldX stdout: {result.stdout}")
                if result.returncode != 0:
                    logger.error(f"FoldX stderr: {result.stderr}")

                # FoldX names mutant output as {pdb_id}_Repair_1.pdb
                mutant_pdb = tmpdir / f"{pdb_id}_Repair_1.pdb"
                if not mutant_pdb.exists():
                    # fallback — search for any _1.pdb
                    candidates = list(tmpdir.glob(f"{pdb_id}_Repair_*.pdb"))
                    if not candidates:
                        raise RuntimeError(
                            f"BuildModel output not found in {tmpdir}. "
                            f"FoldX stdout: {result.stdout[-500:]}"
                        )
                    mutant_pdb = candidates[0]

                with open(mutant_pdb, 'r') as f:
                    mutant_content = f.read()

                ddg = self._parse_ddg(tmpdir / f"Dif_{pdb_id}_Repair.fxout")
                logger.info(f"BuildModel successful: ΔΔG = {ddg}")

                return mutant_content, ddg

            except subprocess.TimeoutExpired:
                raise RuntimeError(f"FoldX timeout after {settings.FOLDX_TIMEOUT_SECONDS}s")
            except Exception as e:
                logger.error(f"BuildModel error: {e}")
                raise
    
    
    def _parse_ddg(self, fxout_file: Path) -> float:
        """Parse ΔΔG from FoldX output"""
        if not fxout_file.exists():
            return 0.0
        
        with open(fxout_file, 'r') as f:
            for line in f:
                if 'total energy difference' in line.lower():
                    parts = line.strip().split()
                    try:
                        return float(parts[-1])
                    except (ValueError, IndexError):
                        pass
        
        return 0.0