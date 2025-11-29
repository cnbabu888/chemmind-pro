import time
import json
import os
from retrosynthesis import retro_engine
from optimization import optimization_engine
from prediction import prediction_engine
from llm_service import llm_service

# Mock the LLM service if no key is present to ensure we test the flow
if not os.environ.get("GEMINI_API_KEY"):
    print("NOTE: GEMINI_API_KEY not found. Testing with MOCK responses.\n")

def test_engine(name, func, input_data):
    print(f"--- Testing {name} ---")
    print(f"Input: {input_data}")
    start_time = time.time()
    try:
        result = func(input_data)
        duration = time.time() - start_time
        print(f"Status: Success")
        print(f"Time: {duration:.4f}s")
        # Print a summary of the result structure
        if "trees" in result:
            print(f"Result: Generated Tree with {len(result['trees'])} root(s)")
        elif "conditions" in result:
            print(f"Result: Suggested {len(result['conditions'])} condition(s)")
            print(f"Top Condition: {result['conditions'][0]['catalyst']} / {result['conditions'][0]['solvent']}")
        elif "product" in result:
            print(f"Result: Predicted Product: {result['product']}")
            print(f"Confidence: {result['confidence']}")
        else:
            print("Result: " + str(result)[:100] + "...")
    except Exception as e:
        print(f"Status: Failed")
        print(f"Error: {e}")
    print("\n")

def main():
    # 1. Retrosynthesis: Ibuprofen
    ibuprofen_smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
    test_engine("Retrosynthesis Engine", retro_engine.predict_route, ibuprofen_smiles)

    # 2. Optimization: Suzuki Coupling
    # Phenylboronic acid + Bromobenzene
    suzuki_reactants = "OB(O)c1ccccc1.Brc1ccccc1" 
    test_engine("Optimization Engine", optimization_engine.suggest_conditions, suzuki_reactants)

    # 3. Forward Prediction: Amide Coupling
    # Acetic acid + Aniline
    amide_reactants = "CC(=O)O.Nc1ccccc1"
    test_engine("Forward Prediction Engine", prediction_engine.predict_product, amide_reactants)

    # 4. Name to Structure (Simulated via LLM Service wrapper if we had it exposed as a function, 
    # but here we'll just check the service status)
    print("--- Testing Name to Structure Service ---")
    if llm_service.model:
        print("LLM Service: Active (Real AI)")
    else:
        print("LLM Service: Inactive (Mock Mode)")
        print("Skipping direct Name-to-Structure test as it requires the API endpoint logic.")

if __name__ == "__main__":
    main()
