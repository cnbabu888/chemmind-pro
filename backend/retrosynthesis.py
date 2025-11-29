import json
import uuid
from typing import Dict, Any
from llm_service import llm_service

class RetrosynthesisEngine:
    def __init__(self):
        pass

    def predict_route(self, target_smiles: str) -> Dict[str, Any]:
        """
        Predicts a synthetic route for the given SMILES using AI.
        """
        if not llm_service.gemini_model and not llm_service.openai_client:
            return self._get_mock_response(target_smiles)

        prompt = f"""
        Act as an expert organic chemist. Perform a retrosynthetic analysis for the molecule with SMILES: "{target_smiles}".
        Provide a plausible synthetic route breaking it down into commercially available precursors.
        
        Return ONLY a JSON object representing the reaction tree. The format must be exactly:
        {{
            "id": "root",
            "smiles": "{target_smiles}",
            "is_chemical": true,
            "children": [
                {{
                    "id": "reaction_1",
                    "is_reaction": true,
                    "is_chemical": false,
                    "metadata": {{
                        "reaction_name": "Name of reaction (e.g. Amide Coupling)",
                        "template_score": 0.95
                    }},
                    "children": [
                        {{
                            "id": "reactant_1",
                            "smiles": "SMILES of reactant 1",
                            "is_chemical": true,
                            "in_stock": true,
                            "children": []
                        }},
                        {{
                            "id": "reactant_2",
                            "smiles": "SMILES of reactant 2",
                            "is_chemical": true,
                            "in_stock": true,
                            "children": []
                        }}
                    ]
                }}
            ]
        }}
        Do not include markdown formatting like ```json. Just the raw JSON string.
        """

        try:
            response_text = llm_service.generate_response(prompt)
            if not response_text:
                return self._get_mock_response(target_smiles)
            
            # Clean up potential markdown code blocks if the model ignores instruction
            cleaned_text = response_text.replace("```json", "").replace("```", "").strip()
            tree_data = json.loads(cleaned_text)
            
            return {
                "target": target_smiles,
                "trees": [tree_data],
                "metadata": {
                    "engine": "Gemini Pro (AIZynthFinder Logic)",
                    "version": "2.0.0-AI",
                    "execution_time": "AI-Generated"
                }
            }
        except Exception as e:
            print(f"Retrosynthesis Error: {e}")
            return self._get_mock_response(target_smiles)

    def _get_mock_response(self, target_smiles: str) -> Dict[str, Any]:
        # Fallback to original mock logic if AI fails or no key
        return {
            "target": target_smiles,
            "trees": [
                {
                    "id": str(uuid.uuid4()),
                    "smiles": target_smiles,
                    "is_chemical": True,
                    "children": [
                        {
                            "id": str(uuid.uuid4()),
                            "is_reaction": True,
                            "is_chemical": False,
                            "metadata": {
                                "reaction_name": "Mock Synthesis (AI Unavailable)",
                                "template_score": 0.5
                            },
                            "children": [
                                {
                                    "id": str(uuid.uuid4()),
                                    "smiles": "C1=CC=CC=C1", # Benzene dummy
                                    "is_chemical": True,
                                    "in_stock": True,
                                    "children": []
                                }
                            ]
                        }
                    ]
                }
            ],
            "metadata": {
                "engine": "Mock Engine (Fallback)",
                "version": "1.0.0",
                "execution_time": "0.0s"
            }
        }

# Global instance
retro_engine = RetrosynthesisEngine()
