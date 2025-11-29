import json
from typing import Dict, Any
from llm_service import llm_service
from rdkit import Chem

class SafetyEngine:
    def __init__(self):
        pass

    def analyze_safety(self, smiles: str) -> Dict[str, Any]:
        """
        Analyzes the safety profile of a molecule using AI.
        Returns toxicity, GHS hazards, and handling precautions.
        """
        alerts = []
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                # Explosive/Reactive patterns
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[N+](=O)[O-]")): alerts.append("Nitro Group (Potential Explosive)")
                if mol.HasSubstructMatch(Chem.MolFromSmarts("N=[N+]=[N-]")): alerts.append("Azide (Explosive)")
                if mol.HasSubstructMatch(Chem.MolFromSmarts("NN")): alerts.append("Hydrazine Derivative (Toxic/Reactive)")
                if mol.HasSubstructMatch(Chem.MolFromSmarts("[O,N]-[O,N]")): alerts.append("Peroxide/N-O bond (Reactive)")
                if mol.HasSubstructMatch(Chem.MolFromSmarts("C#N")): alerts.append("Nitrile (Toxic)")
                if mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)H")): alerts.append("Aldehyde (Reactive)")
        except:
            pass

        if not llm_service.gemini_model and not llm_service.openai_client:
            mock = self._get_mock_response(smiles)
            mock["structural_alerts"] = alerts
            return mock

        prompt = f"""
        Act as an expert chemical safety officer. Analyze the safety profile of the molecule with SMILES: "{smiles}".
        Provide a detailed safety assessment including predicted toxicity, GHS hazards, and handling precautions.
        
        Return ONLY a JSON object. The format must be exactly:
        {{
            "molecule_name": "Common Name",
            "ghs_hazards": [
                {{"code": "H300", "description": "Fatal if swallowed", "category": "Acute Toxicity"}},
                {{"code": "H225", "description": "Highly flammable liquid and vapour", "category": "Flammability"}}
            ],
            "toxicity": {{
                "ld50_rat_oral": "Predicted LD50 value (e.g. 200 mg/kg)",
                "carcinogenicity": "Assessment (e.g. Suspected carcinogen)",
                "mutagenicity": "Assessment"
            }},
            "handling": [
                "Wear protective gloves/protective clothing/eye protection/face protection.",
                "Keep away from heat/sparks/open flames/hot surfaces."
            ],
            "ppe": ["Gloves", "Safety Goggles", "Lab Coat", "Fume Hood"],
            "signal_word": "DANGER" or "WARNING"
        }}
        Do not include markdown formatting like ```json. Just the raw JSON string.
        """

        try:
            response_text = llm_service.generate_response(prompt)
            if not response_text:
                return self._get_mock_response(smiles)
            
            cleaned_text = response_text.replace("```json", "").replace("```", "").strip()
            data = json.loads(cleaned_text)
            data["structural_alerts"] = alerts
            return data

        except Exception as e:
            print(f"Safety Analysis Error: {e}")
            return self._get_mock_response(smiles)

    def _get_mock_response(self, smiles: str) -> Dict[str, Any]:
        # Fallback mock logic
        return {
            "molecule_name": "Unknown Compound",
            "ghs_hazards": [
                {"code": "H302", "description": "Harmful if swallowed", "category": "Acute Toxicity"},
                {"code": "H315", "description": "Causes skin irritation", "category": "Skin Corrosion/Irritation"}
            ],
            "toxicity": {
                "ld50_rat_oral": "500 mg/kg (Estimated)",
                "carcinogenicity": "No data available",
                "mutagenicity": "No data available"
            },
            "handling": [
                "Wash hands thoroughly after handling.",
                "Wear protective gloves and eye protection."
            ],
            "ppe": ["Gloves", "Safety Glasses", "Lab Coat"],
            "signal_word": "WARNING"
        }

# Global instance
safety_engine = SafetyEngine()
