import json
import random
from typing import Dict, Any
from llm_service import llm_service

class ForwardPredictionEngine:
    def __init__(self):
        pass

    def predict_product(self, reactants_smiles: str) -> Dict[str, Any]:
        """
        Predicts the product for the given reactants using AI.
        """
        if not llm_service.gemini_model and not llm_service.openai_client:
            return self._get_mock_response(reactants_smiles)

        prompt = f"""
        Act as an expert organic chemist. Predict the major product for the reaction between these reactants: "{reactants_smiles}".
        Reactants are dot-separated.
        
        Return ONLY a JSON object. The format must be exactly:
        {{
            "product": "SMILES of the major product",
            "confidence": 0.95,
            "reaction_type": "Name of Reaction (e.g. Amide Formation)",
            "explanation": "Brief explanation of the mechanism"
        }}
        Do not include markdown formatting like ```json. Just the raw JSON string.
        """

        try:
            response_text = llm_service.generate_response(prompt)
            if not response_text:
                return self._get_mock_response(reactants_smiles)
            
            cleaned_text = response_text.replace("```json", "").replace("```", "").strip()
            data = json.loads(cleaned_text)
            
            return {
                "reactants": reactants_smiles,
                "product": data.get("product", ""),
                "confidence": data.get("confidence", 0.5),
                "reaction_type": data.get("reaction_type", "Unknown"),
                "metadata": {
                    "engine": "Gemini Pro (IBM RXN Logic)",
                    "explanation": data.get("explanation", ""),
                    "attention_map": [0.5, 0.5] # Placeholder
                }
            }

        except Exception as e:
            print(f"Prediction Error: {e}")
            return self._get_mock_response(reactants_smiles)

    def _get_mock_response(self, reactants_smiles: str) -> Dict[str, Any]:
        # Fallback mock logic
        reactants_list = reactants_smiles.split('.')
        return {
            "reactants": reactants_smiles,
            "product": reactants_list[0] if reactants_list else "",
            "confidence": 0.3,
            "reaction_type": "Low Confidence Prediction (AI Unavailable)",
            "metadata": {
                "engine": "Mock Engine (Fallback)",
                "attention_map": [0.5, 0.5]
            }
        }

# Global instance
prediction_engine = ForwardPredictionEngine()
