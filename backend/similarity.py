from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from typing import List, Dict, Any

class SimilarityEngine:
    def __init__(self):
        # Mock database of common molecules for demonstration
        # In a real app, this would be a vector database or SQL query
        self.database = [
            {"name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"},
            {"name": "Paracetamol", "smiles": "CC(=O)NC1=CC=C(O)C=C1"},
            {"name": "Ibuprofen", "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"},
            {"name": "Caffeine", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"},
            {"name": "Benzene", "smiles": "c1ccccc1"},
            {"name": "Toluene", "smiles": "Cc1ccccc1"},
            {"name": "Phenol", "smiles": "Oc1ccccc1"},
            {"name": "Aniline", "smiles": "Nc1ccccc1"},
            {"name": "Acetic Acid", "smiles": "CC(=O)O"},
            {"name": "Ethanol", "smiles": "CCO"},
            {"name": "Salicylic Acid", "smiles": "OC1=CC=CC=C1C(=O)O"},
            {"name": "Benzoic Acid", "smiles": "O=C(O)c1ccccc1"},
            {"name": "Naproxen", "smiles": "COc1ccc2cc(ccc2c1)C(C)C(=O)O"},
            {"name": "Acetaminophen", "smiles": "CC(=O)Nc1ccc(O)cc1"},
            {"name": "Morphine", "smiles": "CN1CCC23C4C1Cc5c2c(c(O)cc5)OC3C(O)C=C4"},
            {"name": "Codeine", "smiles": "CN1CCC23C4C1Cc5c2c(c(OC)cc5)OC3C(O)C=C4"},
            {"name": "Heroin", "smiles": "CC(=O)OC1C=CC2C3Cc4c(c(OC(C)=O)cc4)OC2C1N(C)CC3"},
            {"name": "Penicillin G", "smiles": "CC1(C(N2C(S1)C(C2=O)NC(=O)Cc3ccccc3)C(=O)O)C"},
            {"name": "Amoxicillin", "smiles": "CC1(C(N2C(S1)C(C2=O)NC(=O)C(c3ccc(O)cc3)N)C(=O)O)C"},
            {"name": "Glucose", "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"},
            {"name": "Fructose", "smiles": "OCC(=O)[C@H](O)[C@@H](O)[C@@H](O)CO"},
            {"name": "Cholesterol", "smiles": "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C"},
            {"name": "Testosterone", "smiles": "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"},
            {"name": "Dopamine", "smiles": "NCCc1ccc(O)c(O)c1"},
            {"name": "Serotonin", "smiles": "NCCc1c[nH]c2ccc(O)cc12"},
            {"name": "Adrenaline", "smiles": "CNC[C@H](O)c1ccc(O)c(O)c1"},
            {"name": "ATP", "smiles": "Nc1ncnc2n(cnc12)[C@@H]3O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]3O"},
            {"name": "Water", "smiles": "O"},
            {"name": "Ammonia", "smiles": "N"},
            {"name": "Methane", "smiles": "C"}
        ]
        
        # Pre-calculate fingerprints for the database
        self.db_fps = []
        for item in self.database:
            mol = Chem.MolFromSmiles(item["smiles"])
            if mol:
                fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
                self.db_fps.append((item, fp))

    def find_similar(self, query_smiles: str, threshold: float = 0.3, limit: int = 10) -> List[Dict[str, Any]]:
        """
        Finds molecules similar to the query SMILES using Tanimoto similarity.
        """
        query_mol = Chem.MolFromSmiles(query_smiles)
        if not query_mol:
            return []

        query_fp = AllChem.GetMorganFingerprintAsBitVect(query_mol, 2, nBits=1024)
        
        results = []
        for item, fp in self.db_fps:
            similarity = DataStructs.TanimotoSimilarity(query_fp, fp)
            if similarity >= threshold:
                # Don't include exact matches (similarity = 1.0) unless requested, 
                # but usually we want "similar" not "same". 
                # Let's include everything >= threshold for now.
                results.append({
                    "name": item["name"],
                    "smiles": item["smiles"],
                    "similarity": round(similarity, 3)
                })
        
        # Sort by similarity descending
        results.sort(key=lambda x: x["similarity"], reverse=True)
        
        return results[:limit]

# Global instance
similarity_engine = SimilarityEngine()
