# app/database.py  ← NEW FILE
import asyncpg
import os

db_pool = None

async def init_db_pool():
    global db_pool
    if db_pool is None:
        db_pool = await asyncpg.create_pool(
            host=os.getenv('DB_HOST', 'localhost'),
            port=int(os.getenv('DB_PORT', 5432)),
            database=os.getenv('DB_NAME', 'alphamissenseAPI'),
            user=os.getenv('DB_USER', 'postgres'),
            password=os.getenv('DB_PASSWORD', '12345'),
            min_size=5,
            max_size=20
        )
    return db_pool