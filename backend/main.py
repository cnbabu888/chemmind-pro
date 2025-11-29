from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from pydantic import BaseModel
from dotenv import load_dotenv
import os

# Load environment variables
load_dotenv()

app = FastAPI(title="ChemMind API", version="0.1.0")

# CORS Configuration
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allow all origins for development
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class MoleculeInput(BaseModel):
    smiles: str

@app.get("/")
def read_root():
    return {"message": "Welcome to ChemMind API"}

@app.get("/health")
def health_check():
    return {"status": "ok"}

@app.post("/api/molecule/info")
def get_molecule_info(data: MoleculeInput):
    mol = Chem.MolFromSmiles(data.smiles)
    if not mol:
        return {"error": "Invalid SMILES"}
    
    return {
        "smiles": data.smiles,
        "num_atoms": mol.GetNumAtoms(),
        "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": Chem.rdMolDescriptors.CalcExactMolWt(mol), # Using Exact Mass as requested for "Exact Mass" field
        "exact_mass": Chem.rdMolDescriptors.CalcExactMolWt(mol),
        "mol_wt": Chem.Descriptors.MolWt(mol), # Average Molecular Weight
        "m_z": Chem.rdMolDescriptors.CalcExactMolWt(mol), # Simplified m/z for M+
        "elemental_analysis": "Calculated from Formula" # Placeholder, or implement logic
    }

@app.get("/api/search")
def search_molecule(query: str):
    """
    Mock search endpoint to resolve names to SMILES.
    """
    # Mock database
    common_molecules = {
        "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "benzene": "c1ccccc1",
        "acetic acid": "CC(=O)O",
        "aniline": "Nc1ccccc1"
    }
    
    query_lower = query.lower()
    if query_lower in common_molecules:
        return {"smiles": common_molecules[query_lower], "name": query}
    
    # If it looks like a SMILES, return it
    if Chem.MolFromSmiles(query):
        return {"smiles": query, "name": "Unknown"}
        
    return {"error": "Molecule not found"}

from fastapi.responses import Response
from rdkit.Chem import Draw

@app.get("/api/molecule/image")
def get_molecule_image(smiles: str):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return Response(content=b"", media_type="image/png", status_code=400)
    
    img = Draw.MolToImage(mol)
    import io
    img_byte_arr = io.BytesIO()
    img.save(img_byte_arr, format='PNG')
    return Response(content=img_byte_arr.getvalue(), media_type="image/png")

from retrosynthesis import retro_engine

@app.post("/api/retrosynthesis")
def run_retrosynthesis(data: MoleculeInput):
    result = retro_engine.predict_route(data.smiles)
    return result

from prediction import prediction_engine

@app.post("/api/prediction")
def run_prediction(data: MoleculeInput):
    # For prediction, we expect 'smiles' to contain dot-separated reactants
    result = prediction_engine.predict_product(data.smiles)
    return result

from chat import chat_engine

class ChatInput(BaseModel):
    message: str
    history: list = []

@app.post("/api/chat")
def chat(data: ChatInput):
    response = chat_engine.get_response(data.message, data.history)
    return {"role": "assistant", "content": response}

from optimization import optimization_engine

from llm_service import llm_service

class NameInput(BaseModel):
    name: str

from safety import safety_engine

class SafetyRequest(BaseModel):
    smiles: str

@app.post("/api/safety")
async def analyze_safety(request: SafetyRequest):
    try:
        result = safety_engine.analyze_safety(request.smiles)
        return result
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/molecule/3d")
async def get_molecule_3d(smiles: str):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            AllChem.MMFFOptimizeMolecule(mol)
            mol_block = Chem.MolToMolBlock(mol)
            return {"mol_block": mol_block}
        return {"error": "Invalid SMILES"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/molecule/name_to_structure")
async def name_to_structure(input_data: NameInput):
    """
    Converts a chemical name to SMILES using AI.
    """
    # Local cache for common molecules to speed up testing and reduce LLM calls
    common_molecules = {
        "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "benzene": "c1ccccc1",
        "acetic acid": "CC(=O)O",
        "aniline": "Nc1ccccc1",
        "water": "O",
        "ethanol": "CCO"
    }
    
    name_lower = input_data.name.lower().strip()
    if name_lower in common_molecules:
        return {"smiles": common_molecules[name_lower]}

    if not llm_service.gemini_model and not llm_service.openai_client:
        # Fallback to simple lookup or error
        return {"error": "AI service unavailable"}

    prompt = f"""
    Act as an expert chemist. Convert the chemical name "{input_data.name}" to its standard SMILES string.
    Return ONLY the SMILES string. No other text.
    """
    
    try:
        # Default to Gemini for now, or could make this configurable
        smiles = llm_service.generate_response(prompt, provider="gemini")
        if smiles:
            return {"smiles": smiles.strip()}
        return {"error": "Could not generate SMILES"}
    except Exception as e:
        return {"error": str(e)}

@app.post("/api/molecule/structure_to_name")
async def structure_to_name(input_data: MoleculeInput):
    """
    Converts a SMILES string to a chemical name using AI.
    """
    if not llm_service.gemini_model and not llm_service.openai_client:
        return {"error": "AI service unavailable"}

    prompt = f"""
    Act as an expert chemist. Convert the SMILES string "{input_data.smiles}" to its common or IUPAC chemical name.
    Return ONLY the name. No other text.
    """
    
    try:
        name = llm_service.generate_response(prompt, provider="gemini")
        if name:
            return {"name": name.strip()}
        return {"error": "Could not generate name"}
    except Exception as e:
        return {"error": str(e)}
@app.post("/api/molecule/structure_to_name")
async def structure_to_name(data: MoleculeInput):
    """
    Converts SMILES to IUPAC name using LLM.
    """
    if not llm_service.gemini_model and not llm_service.openai_client:
        return {"error": "AI service unavailable"}
    
    prompt = f"""
    Give the preferred IUPAC name for the molecule with SMILES: "{data.smiles}".
    Return ONLY the name. No markdown, no explanation.
    """
    try:
        name = llm_service.generate_response(prompt, provider="gemini")
        return {"name": name.strip()}
    except Exception as e:
        return {"error": str(e)}

@app.post("/api/molecule/name_to_structure")
async def name_to_structure(request: NameToStructureRequest):
    """
    Converts chemical name to SMILES using LLM + Cache.
    """
    # 1. Check Cache
    normalized_name = request.name.lower().strip()
    if normalized_name in COMMON_MOLECULES:
        return {"smiles": COMMON_MOLECULES[normalized_name]}

    # 2. Use LLM
    if not llm_service.gemini_model and not llm_service.openai_client:
        return {"error": "AI service unavailable"}
    
    prompt = f"""
    Give the SMILES string for the chemical name: "{request.name}".
    Return ONLY the SMILES string. No markdown, no explanation.
    """
    try:
        smiles = llm_service.generate_response(prompt, provider="gemini")
        # Basic validation: check if it looks like SMILES (not perfect)
        smiles = smiles.strip()
        if " " in smiles or len(smiles) < 1:
             return {"error": "Could not generate name"}
        return {"smiles": smiles}
    except Exception as e:
        return {"error": str(e)}

@app.post("/api/molecule/predict_nmr")
async def predict_nmr(input_data: MoleculeInput):
    """
    Predicts 1H NMR shifts using AI.
    """
    return await _predict_nmr_generic(input_data.smiles, "1H")

@app.post("/api/molecule/predict_c13")
async def predict_c13(input_data: MoleculeInput):
    """
    Predicts 13C NMR shifts using AI.
    """
    return await _predict_nmr_generic(input_data.smiles, "13C")

async def _predict_nmr_generic(smiles: str, nucleus: str):
    if not llm_service.gemini_model and not llm_service.openai_client:
        return {"error": "AI service unavailable"}

    prompt = f"""
    Act as an expert spectroscopist. Predict the approximate {nucleus} NMR spectrum for the molecule with SMILES: "{smiles}".
    Return a JSON object with a list of peaks. Each peak should have:
    - "shift": approximate ppm value (float)
    - "multiplicity": s, d, t, q, m, etc. (for 1H) or just 's' (for 13C decoupled)
    - "count": number of nuclei (integer)
    - "assignment": brief description (e.g., "carbonyl carbon")
    
    Return ONLY the JSON object. No markdown formatting.
    """
    
    try:
        response_text = llm_service.generate_response(prompt, provider="gemini")
        response_text = response_text.replace("```json", "").replace("```", "").strip()
        import json
        data = json.loads(response_text)
        return data
    except Exception as e:
        return {"error": str(e)}

@app.post("/api/molecule/3d")
def get_3d_structure(data: MoleculeInput):
    """
    Generates 3D coordinates for a molecule.
    """
    try:
        mol = Chem.MolFromSmiles(data.smiles)
        if not mol:
            return {"error": "Invalid SMILES"}
        
        mol_h = Chem.AddHs(mol)
        Chem.AllChem.EmbedMolecule(mol_h, randomSeed=42)
        Chem.AllChem.MMFFOptimizeMolecule(mol_h)
        
        # Return MolBlock
        return {"molblock": Chem.MolToMolBlock(mol_h)}
    except Exception as e:
        return {"error": str(e)}

class AnalyzeRequest(BaseModel):
    prompt: str
    provider: str = "gemini"
    deep_search: bool = False

@app.post("/api/analyze")
def analyze_prompt(request: AnalyzeRequest):
    response = llm_service.generate_response(request.prompt, request.provider, request.deep_search)
    return {"response": response}

@app.post("/api/optimization")
def optimize_conditions(data: MoleculeInput):
    # We reuse MoleculeInput but interpret 'smiles' as the reaction mixture/reactants
    result = optimization_engine.suggest_conditions(data.smiles)
    return result

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
