import shutil
import tempfile
from collections import OrderedDict

from dockstream.loggers.ligand_preparation_logger import LigandPreparationLogger
from dockstream.loggers.blank_logger import BlankLogger
from dockstream.utils.dockstream_exceptions import LigandPreparationFailed

from dockstream.utils.execute_external.TautEnum import TautEnumExecutor
from dockstream.utils.enums.taut_enum_enums import TautEnumEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum
from dockstream.core.ligand.ligand import Ligand, get_next_enumeration_number_for_ligand
from dockstream.utils.general_utils import gen_temp_file


class TautEnumSmilePreparator:
    """Class that acts as an interface to the "TautEnum" executable prepare and annotate SMILES."""

    def __init__(self, enumerate_protonation: bool, original_enumeration: bool,
                 add_numbers_to_name: bool, prefix_execution=None, binary_location=None):
        self._TE = TautEnumEnum()
        self._LE = LoggingConfigEnum()
        self._logger = LigandPreparationLogger()
        self._logger_blank = BlankLogger()

        self._enumerate_protonation = enumerate_protonation
        self._original_enumeration = original_enumeration
        self._add_numbers_to_name = add_numbers_to_name
        self._prefix_execution = prefix_execution
        self._binary_location = binary_location

        # check, if backend is available
        self._TautEnum_executor = TautEnumExecutor(prefix_execution=self._prefix_execution,
                                                   binary_location=self._binary_location)
        if not self._TautEnum_executor.is_available():
            raise LigandPreparationFailed("Cannot initialize TautEnum backend - abort.")
        self._logger.log(f"Checked taut_enum backend availability (prefix_execution={prefix_execution}).", self._LE.DEBUG)

    def annotate_tautomers(self, ligands: list) -> list:
        """Method to build all the tautomers for the input SMILES."""

        # 1) generate temporary folder and files
        tmp_output_dir_path = tempfile.mkdtemp()
        tmp_input_smiles_path = gen_temp_file(suffix=".smi", dir=tmp_output_dir_path)
        tmp_output_smiles_path = gen_temp_file(suffix=".smi", dir=tmp_output_dir_path)

        # 2) save the SMILES
        original_smiles = []
        with open(tmp_input_smiles_path, 'w') as f:
            for lig in ligands:
                f.write(lig.get_smile() + ' ' + str(lig.get_ligand_number()) + "\n")
                original_smiles.append(lig.get_original_smile())
        self._logger.log(f"Wrote {len(ligands)} smiles to file {tmp_input_smiles_path} for taut_enum input.", self._LE.DEBUG)

        # 3) run "TautEnum"
        list_args = [self._TE.TAUTENUM_I, tmp_input_smiles_path,
                     self._TE.TAUTENUM_O, tmp_output_smiles_path]
        if self._enumerate_protonation:
            list_args.append(self._TE.TAUTENUM_ENUM_PROTO)
        if self._original_enumeration:
            list_args.append(self._TE.TAUTENUM_ORI_ENUM)
        if self._add_numbers_to_name:
            list_args.append(self._TE.TAUTENUM_ADD_NUMBERS)
        result = self._TautEnum_executor.execute(command=self._TE.TAUTENUM,
                                                 arguments=list_args,
                                                 check=False)
        self._logger.log(f"Executed taut_enum (output file: {tmp_output_smiles_path}).", self._LE.DEBUG)

        # 4) generate a dictionary, where the ligand number is matched to the name
        names_dict = self._get_name_dict(ligands)

        # 5) load and return the smiles; taut_enum output: "COc1cc(c(c(c1OC)OC)Cl)Cc2nc3c(ncnc3n2CCCC#C)N 0_1"
        taut_smiles = []
        taut_identity = []
        with open(tmp_output_smiles_path, 'r') as f:
            for line in f:
                self._logger_blank.log(line.rstrip("\n"), self._LE.DEBUG)
                line = line.strip().split(sep=' ')
                taut_smiles.append(line[0])
                taut_identity.append(line[1])

        buffer = OrderedDict()
        for old_lig in ligands:
            key = str(old_lig.get_ligand_number())
            buffer[key] = {"lig_list": [], "old_lig": old_lig}
        for smile, total_id in zip(taut_smiles, taut_identity):
            total_id_parts = total_id.split('_')
            ligand_number = int(total_id_parts[0])
            for old_id in buffer.keys():
                if old_id == str(ligand_number):
                    matched_list = buffer[old_id]["lig_list"]
                    old_lig = buffer[old_id]["old_lig"]

                    matched_list.append(Ligand(smile=smile,
                                               original_smile=old_lig.get_original_smile(),
                                               ligand_number=ligand_number,
                                               enumeration=len(matched_list),
                                               molecule=None,
                                               mol_type=None,
                                               name=old_lig.get_name()))
        result_list = []
        for key in buffer.keys():
            old_lig = buffer[key]["old_lig"]
            matched_list = buffer[key]["lig_list"]
            if len(matched_list) == 0:
                result_list.append(old_lig)
                continue
            for new_lig in matched_list:
                result_list.append(new_lig)

        # 5) remove temporary files
        shutil.rmtree(tmp_output_dir_path)

        return result_list

    def _get_name_dict(self, ligands: list):
        r_dict = {}
        for lig in ligands:
            r_dict[str(lig.get_ligand_number())] = lig.get_name()
        self._logger.log(f"Using the following dictionary to match the ligand number with the names:\n{r_dict}.", self._LE.DEBUG)
        return r_dict
