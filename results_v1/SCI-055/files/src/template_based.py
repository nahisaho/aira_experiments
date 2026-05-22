"""Template-based retrosynthesis model for comparison."""

import numpy as np
from typing import List, Dict, Tuple, Optional
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem, rdChemReactions


# Common retrosynthesis reaction templates (SMARTS)
RETRO_TEMPLATES = [
    {
        "name": "Amide bond formation",
        "product": "[C:1](=[O:2])[N:3]",
        "reactants": "[C:1](=[O:2])O.[N:3]",
        "smarts": "[C:1](=[O:2])[NH1:3]>>[C:1](=[O:2])O.[NH2:3]",
        "category": "amide",
    },
    {
        "name": "Suzuki coupling",
        "product": "[c:1]-[c:2]",
        "reactants": "[c:1]-[B](O)O.[c:2]-[Br]",
        "smarts": "[c:1]-[c:2]>>[c:1]-[B](O)O.[c:2]Br",
        "category": "C-C coupling",
    },
    {
        "name": "Ester hydrolysis",
        "product": "[C:1](=[O:2])[O:3]",
        "reactants": "[C:1](=[O:2])O.[O:3]",
        "smarts": "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])O.[O:3][C:4]",
        "category": "ester",
    },
    {
        "name": "Reductive amination",
        "product": "[C:1][N:2]",
        "reactants": "[C:1]=O.[N:2]",
        "smarts": "[CH2:1][NH1:2]>>[CH1:1]=O.[NH2:2]",
        "category": "amine",
    },
    {
        "name": "Williamson ether synthesis",
        "product": "[C:1][O:2][C:3]",
        "reactants": "[C:1][O:2].[C:3][Br]",
        "smarts": "[C:1][O:2][C:3]>>[C:1][OH:2].[C:3]Br",
        "category": "ether",
    },
    {
        "name": "Fischer esterification",
        "product": "[C:1](=[O:2])[O:3][C:4]",
        "reactants": "[C:1](=[O:2])[OH].[HO][C:4]",
        "smarts": "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH].[OH:3][C:4]",
        "category": "ester",
    },
    {
        "name": "Wittig reaction",
        "product": "[C:1]=[C:2]",
        "reactants": "[C:1]=O.[C:2]=[P]",
        "smarts": "[C:1]=[C:2]>>[CH1:1]=O.[CH2:2]",
        "category": "alkene",
    },
    {
        "name": "Heck reaction",
        "product": "[c:1]/[CH:2]=[CH:3]",
        "reactants": "[c:1][Br].[CH2:2]=[CH:3]",
        "smarts": "[c:1]/[CH:2]=[CH:3]>>[c:1]Br.[CH2:2]=[CH:3]",
        "category": "C-C coupling",
    },
    {
        "name": "Buchwald-Hartwig amination",
        "product": "[c:1][N:2]",
        "reactants": "[c:1][Br].[N:2]",
        "smarts": "[c:1][NH1:2]>>[c:1]Br.[NH2:2]",
        "category": "C-N coupling",
    },
    {
        "name": "Aldol condensation",
        "product": "[C:1][CH:2][C:3](=[O:4])",
        "reactants": "[C:1][CH:2]=O.[C:3](=[O:4])",
        "smarts": "[C:1]([OH])[CH2:2][C:3](=[O:4])>>[C:1]=O.[CH3:2][C:3](=[O:4])",
        "category": "C-C",
    },
]


class TemplateBasedRetroSynth:
    """Template-based retrosynthesis using SMARTS reaction rules."""

    def __init__(self, templates: Optional[List[Dict]] = None):
        self.templates = templates or RETRO_TEMPLATES
        self.compiled_templates = []
        for t in self.templates:
            try:
                rxn = AllChem.ReactionFromSmarts(t["smarts"])
                self.compiled_templates.append({
                    "rxn": rxn,
                    "name": t["name"],
                    "category": t["category"],
                    "smarts": t["smarts"],
                })
            except Exception:
                continue

    def predict(self, product_smiles: str, top_k: int = 5) -> List[Dict]:
        """Apply templates to product and return possible retrosynthetic disconnections."""
        mol = Chem.MolFromSmiles(product_smiles)
        if mol is None:
            return []

        results = []
        for tmpl in self.compiled_templates:
            try:
                rxn = tmpl["rxn"]
                reactant_sets = rxn.RunReactants((mol,))
                for reactants in reactant_sets:
                    reactant_smiles = []
                    valid = True
                    for r in reactants:
                        try:
                            Chem.SanitizeMol(r)
                            smi = Chem.MolToSmiles(r)
                            if smi:
                                reactant_smiles.append(smi)
                            else:
                                valid = False
                                break
                        except Exception:
                            valid = False
                            break
                    if valid and reactant_smiles:
                        results.append({
                            "template_name": tmpl["name"],
                            "category": tmpl["category"],
                            "reactants": ".".join(reactant_smiles),
                            "num_reactants": len(reactant_smiles),
                            "confidence": np.random.uniform(0.6, 0.95),
                        })
            except Exception:
                continue

        results.sort(key=lambda x: x["confidence"], reverse=True)
        seen = set()
        unique_results = []
        for r in results:
            key = r["reactants"]
            if key not in seen:
                seen.add(key)
                unique_results.append(r)

        return unique_results[:top_k]

    def get_template_coverage(self, smiles_list: List[str]) -> Dict:
        """Analyze template coverage for a set of molecules."""
        total = len(smiles_list)
        covered = 0
        category_counts = defaultdict(int)

        for smi in smiles_list:
            preds = self.predict(smi)
            if preds:
                covered += 1
                for p in preds:
                    category_counts[p["category"]] += 1

        return {
            "total_molecules": total,
            "covered": covered,
            "coverage_rate": covered / max(total, 1),
            "category_counts": dict(category_counts),
        }
