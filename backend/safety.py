import json
from typing import Dict, Any
from llm_service import llm_service

class SafetyEngine:
    def __init__(self):
        pass

    def analyze_safety(self, smiles: str) -> Dict[str, Any]:
        """
        Analyzes the safety profile of a molecule using AI.
        Returns toxicity, GHS hazards, and handling precautions.
        """
        if not llm_service.model:
            return self._get_mock_response(smiles)

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
