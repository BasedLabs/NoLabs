import os
import shutil
import tempfile
import openbabel.openbabel as obab
from rdkit import Chem

from dockstream.utils.dockstream_exceptions import TargetPreparationFailed

from dockstream.core.target_preparator import TargetPreparator
from dockstream.utils.execute_external.rDock import rDockExecutor
from dockstream.utils.enums.rDock_enums import rDockExecutablesEnum, rDockResultKeywordsEnum, rDockRbcavityOutputEnum, rDockTargetPreparationEnum
from dockstream.containers.target_preparation_container import TargetPreparationContainer

from dockstream.utils.general_utils import *


class rDockTargetPreparator(TargetPreparator):
    """Class that deals with all the target preparatory steps needed before docking using "rDock" can commence."""

    def __init__(self, conf: TargetPreparationContainer, target, run_number=0):
        self._TP = rDockTargetPreparationEnum()
        self._EE = rDockExecutablesEnum()
        self._GK = rDockResultKeywordsEnum()
        self._RCO = rDockRbcavityOutputEnum()

        # invoke base class's constructor first
        super().__init__(conf=conf, run_number=run_number)

        # check, whether the backend run specified is an "rDock" one
        if self._run_parameters[self._TP.RUNS_BACKEND] != self._TP.RUNS_BACKEND_RDOCK:
            raise TargetPreparationFailed("Tried to make an rDock preparation with different backend specification.")

        # treat the target: either load a file or store the molecule internally
        if isinstance(target, str):
            if os.path.isfile(target):
                _, file_extension = os.path.splitext(target)
                if file_extension == ".pdb":
                    self._target = Chem.MolFromPDBFile(target, sanitize=False)
                elif file_extension == ".mol2":
                    self._target = Chem.MolFromMol2File(target, sanitize=False)
                else:
                    raise TargetPreparationFailed("Input target file past must end on either \".pdb\" or \".mol2\".")
                self._logger.log(f"Target preparation: File {target} loaded.", self._TL.DEBUG)
            else:
                raise TargetPreparationFailed("Input target file does not exist.")
        elif isinstance(target, Chem.Mol):
            self._target = target
        else:
            raise TargetPreparationFailed("Constructor only accepts a Mol (RDkit) object or a file path.")
        self._logger.log("Stored target as RDkit molecule.", self._TL.DEBUG)

        # initialize the executor for all "rDock" related calls and also check if it is available
        prefix_execution = nested_get(self._run_parameters, [self._TP.RUNS_PARAM,
                                                             self._TP.RUNS_PARAM_PREFIX_EXECUTION],
                                      default=None)
        binary_location = nested_get(self._run_parameters, [self._TP.RUNS_PARAM,
                                                            self._TP.RUNS_PARAM_BINARY_LOCATION],
                                     default=None)
        self._rDock_executor = rDockExecutor(prefix_execution=prefix_execution, binary_location=binary_location)
        if not self._rDock_executor.is_available():
            raise TargetPreparationFailed("Cannot initialize rDock preparator, as rDock backend is not available - abort.")
        self._rDock_executor.set_env_vars()
        self._logger.log(f"Checked rDock backend availability (prefix_execution={prefix_execution}).", self._TL.DEBUG)

    def _export_as_mol2(self, path):
        # rDock expects a Tripos Mol2 file - BUT: there are many different implementations and
        # the RDkit developers decided to go for the "Corina" specification, but for rDock, (modified) Sybyl atom types
        # are required (partial charges are ignored by the scoring functions, though)
        # fun fact: Sybyl atom types are notoriously ill-defined and even differ from tool to tool from the same company
        # --> bottom line: use openbabel instead of RDkit

        # generate temporary copy in Mol2 file
        temp_target_pdb = gen_temp_file(suffix=".pdb")
        Chem.MolToPDBFile(mol=self._target, filename=temp_target_pdb)

        # kekulization of bonds leads to problems with PDB conversion of aromatic side-chains - ignore for now
        # set up conversion; do not forget to add hydrogens (all of them, rDock will remove non-polar ones); translate
        conv = obab.OBConversion()
        conv.SetInAndOutFormats("pdb", "mol2")
        buffer_molecule = obab.OBMol()
        conv.ReadFile(buffer_molecule, temp_target_pdb)
        buffer_molecule.DeleteHydrogens()
        buffer_molecule.AddHydrogens()
        conv.WriteFile(buffer_molecule, path)

        # clean up the temporary file
        if os.path.exists(temp_target_pdb):
            os.remove(temp_target_pdb)
        self._logger.log(f"Exported target as MOL2 file {path}.", self._TL.DEBUG)

    def _parse_rbcavity_output(self, input_str: str) -> dict:
        # assume parameter "input_str" contains text printed to the standard output by the "rbcavity" executable
        result = {}
        lines = [line.strip() for line in input_str.split("\n")]

        # total volume, example: Total volume 1025.75 A^3
        line = [line for line in lines if "Total volume" in line][0]
        result[self._GK.SPECIFYCAVITY_METADATA_TOTALVOLUME] = float(line.split()[2])

        # cavity 1, example: Cavity #1	Size=8206 points; Vol=1025.75 A^3; Min=(-5,1.5,18.5); Max=(13.5,20,31.5);
        # Center=(4.5011,9.94851,24.5739); Extent=(18.5,18.5,13)
        line = [line for line in lines if "Cavity #1" in line][0]
        parts = line.split()
        result[self._GK.SPECIFYCAVITY_METADATA_SIZEINPOINTS] = int(parts[2].split(sep='=')[1])
        self._logger.log(f"Cavity total volume [A^3]: {result[self._GK.SPECIFYCAVITY_METADATA_TOTALVOLUME]}.", self._TL.DEBUG)
        self._logger.log(f"Cavity size in points: {result[self._GK.SPECIFYCAVITY_METADATA_SIZEINPOINTS]}.", self._TL.DEBUG)
        return result

    def _update_prm_file(self, prm_path: str, receptor_mol2_path: str, ref_ligand_sdf_path=None):
        # this function takes a prm file, replaces all occurrences of certain placeholders and writes the updated
        # version to the specified (temporary) path

        with open(prm_path, "r") as prm_file:
            prm_string = prm_file.read()

        # replace the path to the receptor MOL2 file (this is mandatory)
        prm_string = prm_string.replace(self._TP.STRING_RECEPTOR_MOL2_PATH,
                                        receptor_mol2_path)

        # replace the path to the reference ligand SDF file, if required
        if self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD] == self._TP.CAVITY_METHOD_REFERENCE:
            prm_string = prm_string.replace(self._TP.STRING_REFERENCE_LIGAND_SDF_PATH,
                                            ref_ligand_sdf_path)

        with open(prm_path, "wt") as prm_file:
            prm_file.write(prm_string)

    def _cavity_by_reference(self, receptor_mol2_path: str, prm_path: str, folder: str):
        # load the reference file (PDB or SDF)
        ref_format = self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_FORMAT].upper()
        if ref_format == self._TP.CAVITY_REFERENCE_FORMAT_PDB:
            ref_mol = Chem.MolFromPDBFile(self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_PATH],
                                          sanitize=True)
        elif ref_format == self._TP.CAVITY_REFERENCE_FORMAT_SDF:
            mol_supplier = Chem.SDMolSupplier(self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_PATH])
            for mol in mol_supplier:
                ref_mol = mol
        else:
            raise TargetPreparationFailed("Specified format not supported.")

        # write out reference file at temporary location (as SDF, which is required by the "rbcavity" executable)
        ref_ligand_sdf_path = os.path.join(folder, "ref_ligand.sdf")
        sdf_writer = Chem.SDWriter(ref_ligand_sdf_path)
        sdf_writer.write(ref_mol)

        self._update_prm_file(prm_path=prm_path,
                              receptor_mol2_path=receptor_mol2_path,
                              ref_ligand_sdf_path=ref_ligand_sdf_path)

    def specify_cavity(self) -> dict:
        # TODO: implement other cavity specification method(s) than "reference_ligand"
        if self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD] != self._TP.CAVITY_METHOD_REFERENCE:
            raise TargetPreparationFailed("Currently, the rDock backend only supports the reference_ligand cavity specification method.")

        # generate temporary folder
        cur_folder = self._run_parameters[self._TP.RUNS_OUTPUT][self._TP.RUNS_OUTPUT_DIRECTORY]

        # generate a temporary Mol2 file as input for "rbcavity"
        target_mol2_path = os.path.join(cur_folder, "target.mol2")
        self._export_as_mol2(target_mol2_path)

        # copy the specified PRM file into the temporary folder; source can be either user-specified or internal
        prm_input_file = os.path.join(cur_folder, "rbcavity_updated.prm")
        original_prm_file_path = nested_get(self._run_parameters, [self._TP.CAVITY, self._TP.CAVITY_PRMFILE],
                                            default=attach_root_path(self._TP.PRM_DEFAULT_PATH))
        shutil.copyfile(original_prm_file_path, prm_input_file)

        # call the cavity method (only "reference ligand" at the moment)
        self._cavity_by_reference(receptor_mol2_path=target_mol2_path,
                                  prm_path=prm_input_file,
                                  folder=cur_folder)
        self._logger.log(f"Wrote updated PRM file {prm_input_file} - based on template file {original_prm_file_path}.", self._TL.DEBUG)

        # output paths and file names are fixed when using rDock, set them here
        path_cavity_binary = os.path.splitext(prm_input_file)[0] + ".as"
        path_cavity_grid = os.path.splitext(prm_input_file)[0] + "_cav1.grd"

        # set up arguments list and execute
        arguments = [self._EE.RBCAVITY_R, prm_input_file, self._EE.RBCAVITY_D, self._EE.RBCAVITY_WAS]
        result = self._rDock_executor.execute(command=self._EE.RBCAVITY,
                                              arguments=arguments,
                                              check=True)
        self._logger.log(f"Generated cavity binary file {path_cavity_binary} and grid file {path_cavity_grid}.", self._TL.DEBUG)

        if self._RCO.DOCKING_SITE not in result.stdout:
            raise TargetPreparationFailed("".join(["Error occurred when executing rbcavity:\n", result.stdout]))

        # prepare the return dictionary
        dict_return = {self._GK.SPECIFYCAVITY_BINARY_PATH: path_cavity_binary}
        dict_return[self._GK.SPECIFYCAVITY_GRID_PATH] = path_cavity_grid

        # the "rbcavity" executable writes a lot of information to the standard output, such as cavity dimensions
        # -> parse it
        dict_return[self._GK.SPECIFYCAVITY_METADATA] = self._parse_rbcavity_output(input_str=result.stdout)

        # return the paths to the generated files
        return dict_return

    def write_target(self, path):
        pass
