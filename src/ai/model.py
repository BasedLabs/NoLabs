import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer, EsmForProteinFolding
from transformers.models.esm.openfold_utils.protein import to_pdb, Protein as OFProtein
from transformers.models.esm.openfold_utils.feats import atom14_to_atom37

from src.ai.exceptions.model_not_loaded_ex import ModelNotLoadedException


class BaseModel:
    def __init__(self, model_name, gpu):
        self.model_name = model_name
        self.model = None
        self.tokenizer = None
        self.gpu = gpu

    def load_model(self):
        """Load model and tokenizer here"""
        pass

    # Method to get raw model outputs
    def _raw_inference(self, input: str):
        pass

    # Method to return raw outputs in the desired format
    def predict(self, input: str):
        pass


class ClassificationModel(BaseModel):
    def __init__(self, model_name, gpu, task = ""):
        super().__init__(model_name, gpu)

    def set_labels(self, labels):
        self.labels = labels

    def load_model(self):
        self.model = AutoModelForSequenceClassification.from_pretrained(self.model_name)
        self.tokenizer = AutoTokenizer.from_pretrained(self.model_name)
        return super().load_model()

    def _raw_inference(self, sequence: str):
        inputs = self.tokenizer(sequence, return_tensors='pt', padding=True, truncation=True)
        with torch.no_grad():
            outputs = self.model(**inputs)
            probabilities = torch.nn.functional.softmax(outputs.logits, dim=-1).tolist()[0]
        return probabilities

    def predict(self, sequence: str):
        if not self.tokenizer or not self.model:
            raise ModelNotLoadedException()

        probabilities = self._raw_inference(sequence)
        prob_table = list(zip(self.labels, probabilities))
        return prob_table

class Folding(BaseModel):

    def __init__(self, model_name, gpu, task = ""):
        super().__init__(model_name, gpu)

        self.model_name = model_name
        self.gpu = gpu

    def load_model(self):
        self.tokenizer = AutoTokenizer.from_pretrained(model_name)
        self.model = EsmForProteinFolding.from_pretrained(model_name, low_cpu_mem_usage=True)
        if self.gpu:
            self.model = self.model.cuda()

        return super().load_model()

    def convert_outputs_to_pdb(self, outputs):
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
                chain_index=outputs["chain_index"][i] if "chain_index" in outputs else None,
            )
            pdbs.append(to_pdb(pred))
        return pdbs


    def _raw_inference(self, sequence: str):

        tokenized_input = tokenizer([sequence], return_tensors="pt", add_special_tokens=False)['input_ids']
        tokenized_input = tokenized_input.cuda()

        output = None

        with torch.no_grad():
            output = model(tokenized_input)

        return output

    def predict(self, sequence: str) -> List[str]:
        if not self.tokenizer or not self.model:
            raise ModelNotLoadedException()
        output = self._raw_inference(sequence)

        pdbs = self.convert_outputs_to_pdb(output)

        return "".join(pdbs)

    

class FunctionPrediction(BaseModel):
    def __init__(self, model_name, gpu):
        super().__init__(model_name, gpu)

    def predict(self, sequence):
        # TODO: complete function prediction from https://github.com/kexinhuang12345/DeepPurpose
        pass
