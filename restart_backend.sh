#!/bin/bash
echo "----------------------------------------"
echo "  ChemMind Backend Restarter"
echo "----------------------------------------"

# 1. Kill process on port 8000 if it exists
PID=$(lsof -ti:8000)
if [ -n "$PID" ]; then
  echo "Stopping existing server (PID: $PID)..."
  kill -9 $PID
  echo "Server stopped."
else
  echo "No server found running on port 8000."
fi

# 2. Start the server
echo "Starting ChemMind Backend..."
cd backend
# Check if virtual environment exists and activate it if so
if [ -d "../.venv" ]; then
    source ../.venv/bin/activate
fi

# Install dependencies just in case (quietly)
pip install -r requirements.txt > /dev/null 2>&1

# Start Uvicorn
python3 -m uvicorn main:app --reload --host 0.0.0.0 --port 8000
