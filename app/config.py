# app/config.py
from pydantic_settings import BaseSettings
from pathlib import Path

class Settings(BaseSettings):
    # Application
    APP_NAME: str = "DataLens FastAPI Mutation Analysis Server"
    DEBUG: bool = True
    HOST: str = "0.0.0.0"
    PORT: int = 8000

    # Claude API
    CLAUDE_API_KEY: str = ""
    
    # FoldX
    FOLDX_PATH: str = "C://libs//foldx5_1Win64"
    
    # Paths
    CACHE_DIR: Path = Path("./pdb_cache")
    LOG_DIR: Path = Path("./logs")
    TEMP_DIR: Path = Path("./temp")
    
    # Cache settings
    CACHE_ENABLED: bool = True
    CACHE_MAX_AGE_DAYS: int = 30
    
    # RCSB PDB
    RCSB_PDB_URL: str = "https://files.rcsb.org/download"
    
    # Database (for AlphaMissense)
    DATABASE_URL: str = "postgresql://postgres:12345@localhost/alphamissenseAPI"
    
    # Processing limits
    MAX_CONCURRENT_REPAIRS: int = 3
    FOLDX_TIMEOUT_SECONDS: int = 2400

    # Constants
    UNIPROT_BASE_URL: str = "https://rest.uniprot.org/uniprotkb"
    
    class Config:
        env_file = ".env"
        case_sensitive = True

settings = Settings()

# Create directories
settings.CACHE_DIR.mkdir(exist_ok=True)
settings.LOG_DIR.mkdir(exist_ok=True)
settings.TEMP_DIR.mkdir(exist_ok=True)