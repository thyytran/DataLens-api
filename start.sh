#!/bin/bash
# Start UniProt Service for DataLens

echo "=== Starting Datalens Service ==="

# Check Python3
if ! command -v python3 &> /dev/null; then
    echo "Error: python3 not found"
    exit 1
fi

# Install dependencies
echo "Installing dependencies..."
pip3 install fastapi uvicorn requests --break-system-packages --quiet

# Make script executable
chmod +x main.py

# Start service
echo ""
echo "Starting FastAPI service on http://localhost:8000"
echo "Press Ctrl+C to stop"
echo ""

python3 main.py

