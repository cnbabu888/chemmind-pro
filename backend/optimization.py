import json
from typing import Dict, Any
from llm_service import llm_service

class OptimizationEngine:
    def __init__(self):
        pass

    def suggest_conditions(self, reactants_smiles: str) -> Dict[str, Any]:
        """
        Suggests optimal conditions for the given reaction using AI.
        """
        if not llm_service.model:
            return self._get_mock_response(reactants_smiles)

        prompt = f"""
        Act as an expert organic chemist. Suggest optimal reaction conditions for a reaction involving these reactants: "{reactants_smiles}".
        Provide 2-3 distinct methods (e.g., different catalysts, solvents, or conditions).
        Include a "reaction_type" classification.
        
        Return ONLY a JSON object. The format must be exactly:
        {{
            "reaction_type": "Name of Reaction",
            "conditions": [
                {{
                    "rank": 1,
                    "catalyst": "Catalyst Name (e.g. Pd(PPh3)4)",
                    "solvent": "Solvent (e.g. THF)",
                    "base": "Base (optional)",
                    "temperature": "Temperature",
                    "time": "Time",
                    "predicted_yield": "Yield %",
                    "source": "Citation (e.g. J. Org. Chem. 2023, 88, 1234)"
                }}
            ]
        }}
        Do not include markdown formatting like ```json. Just the raw JSON string.
        """

        try:
            response_text = llm_service.generate_response(prompt)
            if not response_text:
                return self._get_mock_response(reactants_smiles)
            
            cleaned_text = response_text.replace("```json", "").replace("```", "").strip()
            data = json.loads(cleaned_text)
            return data

        except Exception as e:
            print(f"Optimization Error: {e}")
            return self._get_mock_response(reactants_smiles)

    def _get_mock_response(self, reactants_smiles: str) -> Dict[str, Any]:
        # Fallback mock logic
        return {
            "reaction_type": "Generic Reaction (AI Unavailable)",
            "conditions": [
                {
                    "rank": 1,
                    "catalyst": "None",
                    "solvent": "Ethanol",
                    "base": "None",
                    "temperature": "Reflux",
                    "time": "24 h",
                    "predicted_yield": "Unknown",
                    "source": "N/A"
                }
            ]
        }

# Global instance
optimization_engine = OptimizationEngine()
