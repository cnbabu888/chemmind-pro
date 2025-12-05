from fastapi import APIRouter, HTTPException
from pydantic import BaseModel
from typing import List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

router = APIRouter()

class ReactionRequest(BaseModel):
    reactants: List[str] # List of SMILES
    reaction_type: Optional[str] = "smart_detect" # or specific name like "suzuki"

class ReactionResponse(BaseModel):
    product_smiles: List[str]
    reaction_name: str
    success: bool
    message: str

# Common Reaction SMARTS Library
REACTION_LIBRARY = {
    "amide_coupling": {
        "name": "Amide Coupling",
        "smarts": "[C:1](=[O:2])[O:3].[N:4]>>[C:1](=[O:2])[N:4]" # Simplified Acid + Amine -> Amide
    },
    "suzuki_coupling": {
        "name": "Suzuki Coupling",
        "smarts": "[c:1][X].[c:2]B(O)O>>[c:1]-[c:2]" # Simplified Aryl Halide + Boronic Acid
    },
    "click_chemistry": {
        "name": "Click Reaction (CuAAC)",
        "smarts": "[C:1]#[C:2].[N:3]=[N:4]=[N:5]>>[C:1]1=[C:2][N:3]([N:4]=[N:5]1)" # Alkyne + Azide -> Triazole
    },
    "esterification": {
        "name": "Fischer Esterification",
        "smarts": "[C:1](=[O:2])[O:3].[O:4]>>[C:1](=[O:2])[O:4]" # Acid + Alcohol -> Ester
    }
}

@router.post("/simulate", response_model=ReactionResponse)
async def simulate_reaction(request: ReactionRequest):
    """
    Simulates a chemical reaction using RDKit SMARTS.
    """
    reactants_mols = [Chem.MolFromSmiles(s) for s in request.reactants]
    if any(m is None for m in reactants_mols):
        raise HTTPException(status_code=400, detail="Invalid SMILES in reactants")

    # 1. Identify Reaction Type (Simple Heuristic or Explicit)
    selected_reaction = None
    
    if request.reaction_type in REACTION_LIBRARY:
        selected_reaction = REACTION_LIBRARY[request.reaction_type]
    else:
        # Auto-detect: Try all and see which one works
        for key, rxn_def in REACTION_LIBRARY.items():
            rxn = AllChem.ReactionFromSmarts(rxn_def["smarts"])
            try:
                products = rxn.RunReactants(tuple(reactants_mols))
                if products:
                    selected_reaction = rxn_def
                    break
            except Exception:
                continue
    
    if not selected_reaction:
        return ReactionResponse(
            product_smiles=[],
            reaction_name="Unknown",
            success=False,
            message="No compatible reaction found for these reactants."
        )

    # 2. Run Reaction
    try:
        rxn = AllChem.ReactionFromSmarts(selected_reaction["smarts"])
        products_set = rxn.RunReactants(tuple(reactants_mols))
        
        if not products_set:
             return ReactionResponse(
                product_smiles=[],
                reaction_name=selected_reaction["name"],
                success=False,
                message="Reaction failed to produce products."
            )
            
        # Flatten products and convert to SMILES
        unique_products = set()
        for product_tuple in products_set:
            for mol in product_tuple:
                try:
                    Chem.SanitizeMol(mol)
                    unique_products.add(Chem.MolToSmiles(mol))
                except:
                    pass
                    
        return ReactionResponse(
            product_smiles=list(unique_products),
            reaction_name=selected_reaction["name"],
            success=True,
            message=f"Successfully simulated {selected_reaction['name']}"
        )

    except Exception as e:
        logger.error(f"Reaction simulation error: {e}")
        raise HTTPException(status_code=500, detail=str(e))
