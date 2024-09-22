from typing import List

import torch
from transformers import AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37
from transformers.models.esm.openfold_utils.protein import Protein as OFProtein
from transformers.models.esm.openfold_utils.protein import to_pdb


class Folding:

    def __init__(self):
        self.model_name = "facebook/esmfold_v1"
        self.gpu = torch.cuda.is_available()

    def load_model(self):
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        self.model = EsmForProteinFolding.from_pretrained(
            self.model_name, low_cpu_mem_usage=True
        )
        if self.gpu:
            self.model = self.model.cuda()
        self.model.eval()

    def convert_outputs_to_pdb(self, outputs) -> List[str]:
        final_atom_positions = atom14_to_atom37(outputs["positions"][-1], outputs)
        outputs = {k: v.to("cpu").numpy() for k, v in outputs.items()}
        final_atom_positions = final_atom_positions.cpu().numpy()
        final_atom_mask = outputs["atom37_atom_exists"]
        pdbs = []
        for i in range(outputs["aatype"].shape[0]):
            aa = outputs["aatype"][i]
            pred_pos = final_atom_positions[i]
            mask = final_atom_mask[i]
            resid = outputs["residue_index"][i] + 1
            pred = OFProtein(
                aatype=aa,
                atom_positions=pred_pos,
                atom_mask=mask,
                residue_index=resid,
                b_factors=outputs["plddt"][i],
                chain_index=(
                    outputs["chain_index"][i] if "chain_index" in outputs else None
                ),
            )
            pdbs.append(to_pdb(pred))
        return pdbs

    def _raw_inference(self, sequence: str):

        tokenized_input = self.tokenizer(
            [sequence], return_tensors="pt", add_special_tokens=False
        )["input_ids"]

        if self.gpu:
            tokenized_input = tokenized_input.cuda()

        with torch.no_grad():
            output = self.model(tokenized_input)

        return output

    def predict(self, sequence: str) -> str:
        output = self._raw_inference(sequence)

        pdbs = self.convert_outputs_to_pdb(output)

        return "".join(pdbs)
