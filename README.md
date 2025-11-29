# ChemMind - AI Chemistry Platform

ChemMind is a next-generation AI-powered chemistry assistant designed to help chemists with retrosynthesis, forward prediction, and condition optimization.

## Features
- **Retrosynthesis**: Visual reaction trees for synthetic route planning.
- **Forward Prediction**: Predict products from reactants with confidence scores.
- **Condition Optimization**: Suggest optimal catalysts, solvents, and conditions.
- **AI Chat**: Conversational interface for chemistry queries.
- **Reports**: Export reaction data to PDF.

   - Frontend: [http://localhost:3000](http://localhost:3000)
   - Backend API: [http://localhost:8000/docs](http://localhost:8000/docs)

### Manual Setup

#### Backend
1. Navigate to `backend`:
   ```bash
   cd backend
   ```
2. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Run server:
   ```bash
   uvicorn main:app --reload
   ```

#### Frontend
1. Navigate to `frontend`:
   ```bash
   cd frontend
   ```
2. Install dependencies:
   ```bash
   npm install
   ```
3. Run dev server:
   ```bash
   npm run dev
   ```

## Tech Stack
- **Frontend**: Next.js, Tailwind CSS, ReactFlow, Lucide Icons
- **Backend**: FastAPI, RDKit, Pydantic
- **Deployment**: Docker, Docker Compose
