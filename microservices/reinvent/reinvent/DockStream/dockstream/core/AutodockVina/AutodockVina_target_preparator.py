from rdkit import Chem

from dockstream.utils.dockstream_exceptions import TargetPreparationFailed

from dockstream.core.target_preparator import TargetPreparator
from dockstream.utils.execute_external.OpenBabel import OpenBabelExecutor
from dockstream.utils.enums.AutodockVina_enums import AutodockTargetPreparationEnum
from dockstream.utils.enums.OpenBabel_enums import OpenBabelExecutablesEnum
from dockstream.containers.target_preparation_container import TargetPreparationContainer

from dockstream.utils.general_utils import *


class AutodockVinaTargetPreparator(TargetPreparator):
    """Class that deals with all the target preparatory steps needed before docking using "Autodock Vina" can commence.

    Note: AutoDockTools recommends to calculate Kollmann (QM-derived but templated for amino acids) charges for the receptor and
    Gasteiger charges for the ligands, but in contrast to Autodock 4, the Vina flavour ignores the charges on the receptor. The
    only thing we need thus to ensure is that (only) the polar hydrogens are present. See FAQs in: http://autodock.scripps.edu/faqs-help/tutorial/using-autodock-with-autodocktools/UsingAutoDockWithADT_v2e.pdf"""

    def __init__(self, conf: TargetPreparationContainer, target, run_number=0):
        self._TP = AutodockTargetPreparationEnum()
        self._EE = OpenBabelExecutablesEnum()

        # invoke base class's constructor first
        super().__init__(conf=conf, run_number=run_number)

        # check, whether the backend run specified is an "rDock" one
        if self._run_parameters[self._TP.RUNS_BACKEND] != self._TP.RUNS_BACKEND_AUTODOCKVINA:
            raise TargetPreparationFailed("Tried to make an AutoDock Vina preparation with different backend specification.")

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

        # initialize the executor for all "OpenBabel" related calls and also check if it is available
        # note, that while there is an "OpenBabel" API (python binding) which we also use, the match to the binary
        # options is not trivial; thus, use command-line here
        self._OpenBabel_executor = OpenBabelExecutor()
        if not self._OpenBabel_executor.is_available():
            raise TargetPreparationFailed("Cannot initialize OpenBabel external library, which should be part of the environment - abort.")
        self._logger.log(f"Checked OpenBabel binary availability.", self._TL.DEBUG)

    def _export_as_pdb2pdbqt(self, path):
        # generate temporary copy
        temp_target_pdb = gen_temp_file(suffix=".pdb")
        Chem.MolToPDBFile(mol=self._target, filename=temp_target_pdb)

        # set target pH value that determines the protein's side-chain states
        if in_keys(self._run_parameters, [self._TP.RUNS_PARAM, self._TP.PH]):
            pH = float(self._run_parameters[self._TP.RUNS_PARAM][self._TP.PH])
        else:
            pH = 7.4
            self._logger.log(f"As a specific pH was not specified, the default pH of {pH} will be used.", self._TL.INFO)

        # Note: In contrast to the ligand preparation, we will not use a tree-based flexibility treatment here - thus,
        #       the option "-xr" is used. Partial charges of the receptor are not used in AutoDock Vina.
        arguments = [temp_target_pdb,
                     self._EE.OBABEL_OUTPUT_FORMAT_PDBQT,
                     "".join([self._EE.OBABEL_O, path]),
                     "".join([self._EE.OBABEL_X, self._EE.OBABEL_X_R]),
                     self._EE.OBABEL_P, pH,
                     self._EE.OBABEL_PARTIALCHARGE, self._EE.OBABEL_PARTIALCHARGE_GASTEIGER]
        self._OpenBabel_executor.execute(command=self._EE.OBABEL,
                                         arguments=arguments,
                                         check=False)

        # clean up the temporary file
        if os.path.exists(temp_target_pdb):
            os.remove(temp_target_pdb)
        self._logger.log(f"Exported target as PDBQT file {path}.", self._TL.DEBUG)

    def _log_extract_box(self):
        x_coords, y_coords, z_coords = self._extract_box()
        if x_coords is not None:
            def dig(value):
                return round(value, ndigits=2)
            self._logger.log(f"Ligand from file {self._run_parameters[self._TP.RUNS_PARAM][self._TP.EXTRACT_BOX][self._TP.EXTRACT_BOX_REFERENCE_LIGAND_PATH]} has the following dimensions:",
                             self._TL.INFO)
            self._logger_blank.log(f"X coordinates: min={dig(min(x_coords))}, max={dig(max(x_coords))}, mean={dig(sum(x_coords)/len(x_coords))}",
                                   self._TL.INFO)
            self._logger_blank.log(f"Y coordinates: min={dig(min(y_coords))}, max={dig(max(y_coords))}, mean={dig(sum(y_coords)/len(y_coords))}",
                                   self._TL.INFO)
            self._logger_blank.log(f"Z coordinates: min={dig(min(z_coords))}, max={dig(max(z_coords))}, mean={dig(sum(z_coords)/len(z_coords))}",
                                   self._TL.INFO)

    def _extract_box(self):
        # extracts box suggestions from a reference ligand, which can be added to a AutoDock Vina run
        if in_keys(self._run_parameters, [self._TP.RUNS_PARAM, self._TP.EXTRACT_BOX]):
            if in_keys(self._run_parameters, [self._TP.RUNS_PARAM, self._TP.EXTRACT_BOX, self._TP.EXTRACT_BOX_REFERENCE_LIGAND_PATH]) and \
               in_keys(self._run_parameters, [self._TP.RUNS_PARAM, self._TP.EXTRACT_BOX, self._TP.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT]):

                # load the reference file (PDB or SDF)
                ref_format = self._run_parameters[self._TP.RUNS_PARAM][self._TP.EXTRACT_BOX][self._TP.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT].upper()
                if ref_format == self._TP.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_PDB:
                    ref_mol = Chem.MolFromPDBFile(self._run_parameters[self._TP.RUNS_PARAM][self._TP.EXTRACT_BOX][self._TP.EXTRACT_BOX_REFERENCE_LIGAND_PATH],
                                                  sanitize=True)
                elif ref_format == self._TP.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT_SDF:
                    mol_supplier = Chem.SDMolSupplier(self._run_parameters[self._TP.RUNS_PARAM][self._TP.EXTRACT_BOX][self._TP.EXTRACT_BOX_REFERENCE_LIGAND_PATH])
                    for mol in mol_supplier:
                        ref_mol = mol
                else:
                    raise TargetPreparationFailed("Specified format not supported.")

                # extract coordinates
                x_coords = [atom[0] for atom in ref_mol.GetConformer(0).GetPositions()]
                y_coords = [atom[1] for atom in ref_mol.GetConformer(0).GetPositions()]
                z_coords = [atom[2] for atom in ref_mol.GetConformer(0).GetPositions()]
                return x_coords, y_coords, z_coords
            else:
                self._logger.log(f"In order extract the box, both {self._TP.EXTRACT_BOX_REFERENCE_LIGAND_PATH} and {self._TP.EXTRACT_BOX_REFERENCE_LIGAND_FORMAT} must be defined.")
                return None, None, None

    def specify_cavity(self):
        # write out the input PDB as PDBQT file
        self._export_as_pdb2pdbqt(self._run_parameters[self._TP.RUNS_OUTPUT][self._TP.RECEPTOR_PATH])

        # if there is a reference ligand provided, calculate mean, minimum and maximum coordinates and log out
        self._log_extract_box()

    def write_target(self, path):
        # TODO: move writing functionality here (and for rDock) to this method, respectively
        pass
