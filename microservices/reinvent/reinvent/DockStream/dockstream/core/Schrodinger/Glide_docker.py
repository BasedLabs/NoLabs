import tempfile
import multiprocessing
import os
import time
import gzip
import shutil
from enum import Enum
from typing import Optional, Dict, List, Union, Iterable
from typing_extensions import Literal  # Required for Python 3.7. From 3.8 Literal is in typing.

from copy import deepcopy

import rdkit.Chem as Chem
from pydantic import PrivateAttr, BaseModel, Field

from dockstream.core.docker import Docker, _LE
from dockstream.core.Schrodinger.license_token_guard import SchrodingerLicenseTokenGuard
from dockstream.core.Schrodinger.Glide_result_parser import GlideResultParser
from dockstream.utils.execute_external.Schrodinger import SchrodingerExecutor
from dockstream.utils.enums.ligand_preparation_enum import LigandPreparationEnum
from dockstream.utils.enums.Schrodinger_enums import SchrodingerExecutablesEnum, \
                                                 SchrodingerDockingConfigurationEnum, \
                                                 SchrodingerOutputEnum
from dockstream.utils.general_utils import gen_temp_file

from dockstream.utils.translations.molecule_translator import MoleculeTranslator
from dockstream.utils.files_paths import any_in_file
from dockstream.utils.dockstream_exceptions import DockingRunFailed

_LP = LigandPreparationEnum()
_CE = SchrodingerDockingConfigurationEnum()
_EE = SchrodingerExecutablesEnum()
_ROE = SchrodingerOutputEnum()


class Parallelization(BaseModel):
    number_cores: Optional[int] = Field(default=4)
    max_compounds_per_subjob: Optional[int] = Field(default=0, ge=0)


class AmideMode(str, Enum):
    PENAL = "penal"
    FIXED = "fixed"
    FREE = "free"
    TRANS = "trans"
    GENERALIZED = "generalized"


class GlidePoseOuttype(str, Enum):
    POSEVIEWER = "poseviewer"
    LIGANDLIB = "ligandlib"
    POSEVIEWER_SD = "poseviewer_sd"
    LIGANDLIB_SD = "ligandlib_sd"  # the only supported type for now
    PHASE_SUBSET = "phase_subset"


class GlidePrecision(str, Enum):
    """Glide Precision.

    HTVS (high-throughput virtual screening),
    SP (standard precision),
    XP (extra precision).
    """

    # HTVS and SP docking use the same scoring function,
    # but HTVS reduces the number of intermediate conformations throughout the docking funnel,
    # and also reduces the thoroughness of the final torsional refinement and sampling.
    # The docking algorithm itself is essentially the same.
    #
    # XP does more extensive sampling than SP:
    # it starts with SP sampling before beginning its own anchor-and-grow procedure.
    # XP also employs a more sophisticated scoring function that is "harder" than the SP GlideScore,
    # with greater requirements for ligand-receptor shape complementarity.
    # This weeds out false positives that SP lets through.
    # Because XP can penalize ligands that don't fit well to the particular receptor conformation used,
    # we recommend docking to multiple receptor conformations, if possible.

    SP = "SP"
    HTVS = "HTVS"
    XP = "XP"


class GlideKeywords(BaseModel):

    # These keywords override whatever is in advanced_glide_keywords.maestro_file!
    # Since each of them has a default value, they always override.
    AMIDE_MODE: AmideMode = Field(default=AmideMode.TRANS, title="AMIDE_MODE", description="""Amide bond sampling mode/rotation behavior: "penal" - penalize nonplanar conformation, "free" - vary conformation, "fixed" - retain original conformation, "trans" - allow trans conformation only, "generalized" - use generalized torsion controls defined in torcontrol.txt.""")
    EXPANDED_SAMPLING: bool = Field(default=False, title="EXPANDED_SAMPLING", description="Bypass elimination of poses in rough scoring stage (useful for fragment docking)")
    GRIDFILE: Union[str, List[str]] = Field(title="GRIDFILE", description="Path to grid (.grd or .zip) file")
    LIGANDFILE: str = None  #: Glide docking ligands file name
    NENHANCED_SAMPLING: int = Field(default=1, ge=1, le=4, title="NENHANCED_SAMPLING", description="Expand size of the Glide funnel by N times to process poses from N confgen runs with minor perturbations to the input ligand coordinates")
    POSE_OUTTYPE: GlidePoseOuttype = Field(default=GlidePoseOuttype.LIGANDLIB_SD, title="POSE_OUTTYPE", description='format for file containing docked poses: "poseviewer" for _pv.mae output; "ligandlib" for _lib.mae; similarly "poseviewer_sd" and "ligandlib_sd" for sdf output; "phase_subset" for bypassing _lib or _pv in favor of a Phase subset file.')
    POSES_PER_LIG: int = Field(default=1, ge=1, title="POSES_PER_LIG", description="maximum number of poses to report per each input ligand")
    POSTDOCK_NPOSE: int = Field(default=1, ge=0, title="POSTDOCK_NPOSE", description="maximum number of best-by-Emodel poses to submit to post-docking minimization")
    POSTDOCKSTRAIN: bool = Field(default=False, title="POSTDOCKSTRAIN", description="include strain correction in post-docking score")
    PRECISION: GlidePrecision = Field(default=GlidePrecision.SP, title="PRECISION", description='Glide docking precision: "HTVS" - high-throughput virtual screening, "SP" - standard precision, "XP" - extra precision.')
    REWARD_INTRA_HBONDS: bool = Field(default=False, title="REWARD_INTRA_HBONDS", description="Reward formation of intramolecular hydrogen bonds in the ligand")

    class Config:

        use_enum_values = True
        extra = "allow"


class AdvancedGlideKeywords(BaseModel):
    """Additional Glide keywords, constraints and features can be provided as Maestro file.
    Keywords in Maestro file will have lower priority than GUI keywords,
    and will be overwritten by keywords specified in the GUI.
    If you specify "USE_REF_LIGAND   True" or REF_LIGAND_FILE  keywords in Maestro,
    you should also upload REF_LIGAND_FILE.
    REF_LIGAND_FILE keyword in Maestro file will be replaced by the uploaded file.
    """

    maestro_file: Optional[str] = Field(default=None, title="Maestro .in file", description="Maestro file with additional keywords, constraints and features.")
    REF_LIGAND_FILE: Optional[str] = Field(default=None, title="Reference ligand file", description="Reference ligand file (required if USE_REF_LIGAND is set to True).")


class GlideParameters(BaseModel):
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None
    time_limit_per_compound: Optional[int] = None
    token_guard: Optional[Dict] = None
    parallelization: Optional[Parallelization]
    glide_flags: Optional[Dict]
    glide_keywords: Optional[GlideKeywords]
    advanced_glide_keywords: Optional[AdvancedGlideKeywords]


def stringify(obj):
    """Converts all objects in a dict to str, recursively."""
    if isinstance(obj, dict):
        return {str(key): stringify(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [stringify(value) for value in obj]
    else:
        return str(obj)


def parse_maestro(lines: Iterable[str]) -> Dict[str, Union[str, Dict[str, str]]]:
    """Parses Maestro input, and returns DockStream keywords dict for it."""

    separator3 = "   "
    indent4 = "    "
    block_starters = {
        "[CONSTRAINT_GROUP",
        "[FEATURE",
    }

    # All Glide keywords. Get all keywords with:
    #   $ module load schrodinger
    #   $ glide -k | cut -d' ' -f1 | sed 's/.*/"&"/' | paste -sd , -
    # List keywords, get first word, wrap in quotes, join lines.
    # See:
    #   - https://stackoverflow.com/a/19145499
    #   - https://unix.stackexchange.com/a/251362
    allowed_keywords = {"AMIDE_MODE","AMIDE_TRANS_ALL","AMIDE_TRANSTOL","ASL_RES_INTERACTION","ASLSTRINGS","CALC_INPUT_RMS","CANONICALIZE","COMPRESS_POSES","CORE_ATOMS","CORE_DEFINITION","CORE_FILTER","CORE_POS_MAX_RMSD","CORE_RESTRAIN","CORE_RESTRAIN_V","CORE_SMARTS","CORE_SNAP","CORECONS_FALLBACK","CSV_PROPS_FILE","CV_CUTOFF","DIELMOD","DOCKING_METHOD","DOINTRA","DOINTRA_SCALE","DSCORE_CUTOFF","EPIK_PENALTIES","EXPANDED_SAMPLING","FITDEN","FLEXASL","FORCEFIELD","FORCEPLANAR","GLIDE_CONFGEN_BADDIST2","GLIDE_CONFGEN_EFCUT","GLIDE_CONS_FEAT_FILE","GLIDE_CONS_FINALONLY","GLIDE_CONS_RMETCOORD","GLIDE_CONS_RNOEMAX","GLIDE_CONS_RNOEMIN","GLIDE_CONS_RPOS","GLIDE_CONS_XMETCOORD","GLIDE_CONS_XNOE","GLIDE_CONS_XPOS","GLIDE_CONS_YMETCOORD","GLIDE_CONS_YNOE","GLIDE_CONS_YPOS","GLIDE_CONS_ZMETCOORD","GLIDE_CONS_ZNOE","GLIDE_CONS_ZPOS","GLIDE_DIELCO","GLIDE_ELEMENTS","GLIDE_EXVOL_PENAL_NUM","GLIDE_EXVOL_PENAL_STRENGTH","GLIDE_NTOTALCONS","GLIDE_NUMEXVOL","GLIDE_NUMMETCOORDCONS","GLIDE_NUMMETCOORDSITES","GLIDE_NUMNOECONS","GLIDE_NUMPOSITCONS","GLIDE_NUMUSEXVOL","GLIDE_OUTPUT_USEHTOR","GLIDE_RECEP_ASLSCALE","GLIDE_RECEP_MAESCALE","GLIDE_REFLIG_FORMAT","GLIDE_REXVOL","GLIDE_REXVOLIN","GLIDE_TORCONS_ALLBONDS","GLIDE_TORCONS_IATOMS","GLIDE_TORCONS_JATOMS","GLIDE_TORCONS_KATOMS","GLIDE_TORCONS_LATOMS","GLIDE_TORCONS_PATTERN_INDEX","GLIDE_TORCONS_PATTERNS","GLIDE_TORCONS_SETVAL","GLIDE_TORCONS_VALUES","GLIDE_TORCONSFILE","GLIDE_XEXVOL","GLIDE_XP_NMAXCORE","GLIDE_XP_RMSCUT","GLIDE_YEXVOL","GLIDE_ZEXVOL","GLIDECONS","GLIDECONSATOMS","GLIDECONSFEATATOMS","GLIDECONSFEATHASINCLUDE","GLIDECONSFEATINCLUDE","GLIDECONSFEATINDEX","GLIDECONSFEATPATTERNS","GLIDECONSGROUPNREQUIRED","GLIDECONSNAMES","GLIDECONSUSEMET","GLIDECONSUSESYMATOMS","GLIDERECEPTORSCALECHARGES","GLIDERECEPTORSCALERADII","GLIDESCORUSEMET","GLIDEUSEALLEXVOL","GLIDEUSECONSFEAT","GLIDEUSECONSFEATINDEX","GLIDEUSECONSGROUPINDEX","GLIDEUSECONSLABELS","GLIDEUSEXVOL","GLIDEUSEXVOLNAMES","GLIDEXVOLNAMES","GRID_CENTER","GRID_CENTER_ASL","GRIDFILE","GSCORE","GSCORE_CUTOFF","HAVEGLIDECONSFEAT","HBOND_ACCEP_HALO","HBOND_CONSTRAINTS","HBOND_CUTOFF","HBOND_DONOR_AROMH","HBOND_DONOR_AROMH_CHARGE","HBOND_DONOR_HALO","INCLUDE_INPUT_CONF","INCLUDE_INPUT_RINGS","INNERBOX","JOBNAME","KEEP_SUBJOB_POSES","KEEPRAW","KEEPSKIPPED","LIG_CCUT","LIG_MAECHARGES","LIG_VSCALE","LIGAND_END","LIGAND_INDEX","LIGAND_MOLECULE","LIGAND_START","LIGANDFILE","LIGANDFILES","LIGFORMAT","LIGPREP","LIGPREP_ARGS","MACROCYCLE","MACROCYCLE_OPTIONS","MAX_ITERATIONS","MAXATOMS","MAXKEEP","MAXREF","MAXROTBONDS","METAL_CONSTRAINTS","METAL_CUTOFF","METCOORD_CONSTRAINTS","METCOORD_SITES","NENHANCED_SAMPLING","NMAXRMSSYM","NOE_CONSTRAINTS","NOSORT","NREPORT","NREQUIRED_CONS","OUTERBOX","OUTPUTDIR","PAIRDISTANCES","PEPTIDE","PHASE_DB","PHASE_NCONFS","PHASE_SUBSET","POSE_DISPLACEMENT","POSE_HTORSION","POSE_OUTTYPE","POSE_RMSD","POSES_PER_LIG","POSIT_CONSTRAINTS","POSTDOCK","POSTDOCK_ITMAX","POSTDOCK_NPOSE","POSTDOCK_SCITMAX","POSTDOCK_XP_DELE","POSTDOCKCG","POSTDOCKLIGMIN","POSTDOCKSTRAIN","PRECISION","PREMIN","PREMINCG","PREMINELEC","PREMINITMAX","RADIUS_RES_INTERACTION","REC_MAECHARGES","RECEP_CCUT","RECEP_FILE","RECEP_VSCALE","REF_LIGAND_FILE","REFINDEX","REPORT_CPU_TIME","REWARD_INTRA_HBONDS","RINGCONFCUT","RINGONFLY","SAMPLE_N_INVERSIONS","SAMPLE_RINGS","SCORE_INPUT_POSE","SCORE_MINIMIZED_INPUT_POSE","SCORING_CUTOFF","SHAPE_ATOMS","SHAPE_RESTRAIN","SHAPE_TYPING","SKIP_EPIK_METAL_ONLY","STRAIN_GSFACTOR","STRAIN_GSTHRESH","STRAINELEC","SUBSTRATE_PENAL_FILE","USE_CONS","USE_REF_LIGAND","USECOMPMAE","USEFLEXASL","USEFLEXMAE","WRITE_CSV","WRITE_RES_INTERACTION","WRITE_TIMINGS_CSV","WRITE_XP_DESC","WRITEREPT"}

    result = {}
    current_block = None
    for linenum, line in enumerate(lines):
        if any(line.startswith(starter) for starter in block_starters):
            # Block start.
            current_block = line.strip()
            result[current_block] = {}
        elif line.strip() == "":
            # Empty line: close current block if any is open, and skip the line.
            current_block = None
        elif line.startswith(indent4):
            # Indented line inside the block.
            if current_block is None:
                raise ValueError(f"Unexpected indent outside of block for line {linenum}: {line}")
            kw, value = line.strip().split(sep=separator3, maxsplit=1)
            result[current_block][kw] = value.strip("\"")
        elif any(line.startswith(kw) for kw in allowed_keywords):
            # Ordinary keywords.
            kw, value = line.strip().split(sep=separator3, maxsplit=1)
            result[kw] = value
        else:
            raise ValueError(f"Unexpected line {linenum}: {line}")

    return result


class Glide(Docker, BaseModel):
    """Interface to the "Glide" backend."""

    backend: Literal["Glide"] = "Glide"
    parameters: GlideParameters = GlideParameters()

    class Config:
        underscore_attrs_are_private = True

    _Schrodinger_executor: SchrodingerExecutor = PrivateAttr()
    _token_guard: SchrodingerLicenseTokenGuard = PrivateAttr(default=None)
    _execution_result = PrivateAttr()  # "Glide" specific return stuff.

    def __init__(self, **data):
        super().__init__(**data)

        self._Schrodinger_executor = SchrodingerExecutor(prefix_execution=self.parameters.prefix_execution,
                                                         binary_location=self.parameters.binary_location)
        if not self._Schrodinger_executor.is_available():
            raise DockingRunFailed("Cannot initialize Glide docker, as Schrodinger backend is not available - abort.")
        self._logger.log(f"Checked Glide backend availability (prefix_execution={self.parameters.prefix_execution}, binary_location={self.parameters.binary_location}).",
                         _LE.DEBUG)

        # if necessary parameters are specified, initialize the token guard
        if self.parameters.token_guard is not None:
            self._token_guard = SchrodingerLicenseTokenGuard(**self.parameters.token_guard)

        self._execution_result = None  # "Glide" specific return stuff.

    def _apply_token_guard(self):
        if self._token_guard is not None:
            self._token_guard.guard()

    def _get_score_from_conformer(self, conformer):
        return float(conformer.GetProp(_ROE.GLIDE_DOCKING_SCORE))

    def add_molecules(self, molecules: list):
        """This method overrides the parent class, docker.py add_molecules method. This method appends prepared
        ligands to a list for subsequent docking. Note, that while internally we will store the ligands for "Glide"
        in RDkit format, they will need to be written out in Schrodinger's MAE format for subsequent docking

        :param molecules: A list that is to contain all prepared ligands for subsequent docking
        :type molecules: list
        :raises NotImplementedError: Each backend must override the parent class, docker.py add_molecules method.
            Inability to do so or a bug causing incorrect implementation will raise a NotImplementedError
        """
        # note, that while internally we will store the ligands for "Glide" in RDkit format, they will need to be
        # written out in MAE before docking can commence later
        mol_trans = MoleculeTranslator(self.ligands, force_mol_type=_LP.TYPE_RDKIT)
        mol_trans.add_molecules(molecules)
        self.ligands = mol_trans.get_as_rdkit()
        self._docking_performed = False

    def get_execution_result(self):
        """This method returns execution_result which contains all the information related to the Glide docking run,
        including LigPrep, docking parameters, etc.

        :raises ValueError: command which is a keyword in execution_result determines which docking step to run
            (ex. LigPrep). The command must be a valid term in the internal Schrodinger dictionary of executables
        :return execution_result: Contains strings and lists related to Glide docking configurations
        """
        return self._execution_result

    def _translate_SDF_to_MAE(self, sdf_path, mae_path):
        """As "Glide" is only able to read MAE (Maestro) files, write the ligands out in that format."""
        # call "sdconvert" from Schrodinger's software
        arguments = ["".join([_EE.SDCONVERT_I, _EE.SDCONVERT_FORMAT_SD]), sdf_path,
                     "".join([_EE.SDCONVERT_O, _EE.SDCONVERT_FORMAT_MAE]),mae_path]
        execution_result = self._Schrodinger_executor.execute(command=_EE.SDCONVERT,
                                                              arguments=arguments,
                                                              check=True)

    def _all_keywords(self):
        """Returns joined keywords from JSON and from "advanced" Maestro .in file."""

        keywords = {}

        # Add keywords from maestro file.
        if self.parameters.advanced_glide_keywords is not None and \
                self.parameters.advanced_glide_keywords.maestro_file is not None:
            with open(self.parameters.advanced_glide_keywords.maestro_file, "rt") as f:
                keywords_from_file = parse_maestro(f.readlines())
                keywords.update(keywords_from_file)

        # Add keywords from advanced_glide_keywords
        # (they are keywords with file paths),
        # skipping keywords that are None.
        # Also skip maestro file - that's not a keyword.
        if self.parameters.advanced_glide_keywords is not None:
            adv_kw = stringify({
                k: v
                for k, v in self.parameters.advanced_glide_keywords.dict().items()
                if v is not None and k not in {'maestro_file'}
            })
            keywords.update(adv_kw)

        # Add "ordinary" keywords, overwriting existing/advanced ones.
        json_keywords = stringify(deepcopy(self.parameters.glide_keywords.dict()))
        keywords.update(json_keywords)  # Overwrites any keywords that are already present.
        return keywords

    def _keywords_Maestro_reformat(self, keywords: dict):
        # rewrite keyword input file in Maestro format
        maestro_indent = "    "
        maestro_spacing = "   "

        element_lines = []
        block_lines = []

        for key in keywords.keys():
            if isinstance(keywords[key], str):
                # keyword holds one dictionary (string) only
                element_lines.append(maestro_spacing.join([key, keywords[key] + "\n"]))
            elif isinstance(keywords[key], dict):
                # keyword holds a composite block and has no dictionary (e.g. constraints); note, that these must
                # always be at the end of the file
                block_lines.append("\n" + key + "\n")
                block = keywords[key]
                for key_idx, block_key in enumerate(block.keys()):
                    block_value = block[block_key]

                    # if this is a value in certain blocks, put it into double quotation marks as spaces are present
                    if any([x in key for x in _CE.GLIDE_INPUTBLOCK_VALUEQUOTED]):
                        block_value = '"' + block_value + '"'
                    line = maestro_indent + maestro_spacing.join([block_key, block_value])

                        # add comma to block definition, if there are more lines to come and the block requires it
                        # note, that not all blocks in GLIDE require this; in some cases, the comma is already part of
                        # the line (then skip it!)
                    if any([x in key for x in _CE.GLIDE_INPUTBLOCK_COMMASEPARATED]):
                        if (key_idx + 1) < len(block) and line[-1] != ',':
                            line = line + ','

                    block_lines.append(line + "\n")
            else:
                raise Exception(f"Cannot handle type {type(keywords[key])} in keyword file specification, only use strings and blocks.")

        return element_lines, block_lines

    def _write_keywords_to_file(self, keywords: dict, path=None) -> str:
        """Function to generate a keyword input file in Maestro format."""
        # TODO: support "ensemble docking" - at the moment only the first gridfile is used
        keywords = deepcopy(keywords)
        gridfiles = keywords[_EE.GLIDE_GRIDFILE]
        if not isinstance(gridfiles, list):
            gridfiles = [gridfiles]
        keywords[_EE.GLIDE_GRIDFILE] = gridfiles[0]

        # enforce, that the write-out mode is the one supported
        keywords[_EE.GLIDE_POSE_OUTTYPE] = _EE.GLIDE_POSE_OUTTYPE_LIGANDLIB

        # call a function that returns the input keywords in Maestro format
        element_lines, block_lines = self._keywords_Maestro_reformat(keywords)

        # arrange the elements and blocks
        if path is None:
            path = gen_temp_file(suffix=".in")
        with open(path, mode='w') as f:
            self._logger.log(f"Writing GLIDE input file {path}:\n", _LE.DEBUG)
            for line in element_lines:
                f.write(line)
                self._logger_blank.log(line.rstrip("\n"), _LE.DEBUG)
            for line in block_lines:
                f.write(line)
                self._logger_blank.log(line.rstrip("\n"), _LE.DEBUG)
            self._logger_blank.log("", _LE.DEBUG)
            self._logger.log("--- End file", _LE.DEBUG)
        return path

    def _generate_temporary_input_output_files(self, start_indices, sublists):
        # in case singletons are handed over, wrap them in a list for "zipping" later
        if not isinstance(start_indices, list):
            start_indices = [start_indices]
        if not isinstance(sublists, list):
            sublists = [sublists]

        tmp_output_dirs = []
        tmp_input_mae_paths = []
        tmp_output_sdf_paths = []
        for start_index, sublist in zip(start_indices, sublists):
            # generate temporary input files and output directory
            cur_tmp_output_dir = tempfile.mkdtemp()
            cur_tmp_sdf = gen_temp_file(prefix=str(start_index), suffix=".sdf", dir=cur_tmp_output_dir)
            cur_tmp_mae = gen_temp_file(prefix=str(start_index), suffix=".mae", dir=cur_tmp_output_dir)

            # write-out the temporary input file
            writer = Chem.SDWriter(cur_tmp_sdf)
            one_written = False
            for ligand in sublist:
                # initialize all ligands (as they could have failed)
                if ligand.get_molecule() is not None:
                    mol = deepcopy(ligand.get_molecule())
                    mol.SetProp("_Name", ligand.get_identifier())
                    one_written = True
                    writer.write(mol)
            writer.close()
            if one_written is False:
                if os.path.isdir(cur_tmp_output_dir):
                    shutil.rmtree(cur_tmp_output_dir)
                continue

            # translate the SDF into a MAE file
            self._translate_SDF_to_MAE(sdf_path=cur_tmp_sdf, mae_path=cur_tmp_mae)

            # add the path to which "_dock_subjob()" will write the result SDF
            output_sdf_path = gen_temp_file(prefix=str(start_index), suffix="_result.sdf", dir=cur_tmp_output_dir)
            tmp_output_sdf_paths.append(output_sdf_path)
            tmp_input_mae_paths.append(cur_tmp_mae)
            tmp_output_dirs.append(cur_tmp_output_dir)
        return tmp_output_dirs, tmp_input_mae_paths, tmp_output_sdf_paths

    def _dock(self, number_cores: int):
        # partition ligands into sublists and distribute to processor cores for docking
        start_indices, sublists = self.get_sublists_for_docking(number_cores=number_cores)
        number_sublists = len(sublists)
        number_ligands_per_sublist = len(sublists[0])
        self._logger.log(f"Split ligands into {number_sublists} sublists for docking.", _LE.DEBUG)
        sublists_submitted = 0
        slices_per_iteration = min(number_cores, number_sublists)

        while sublists_submitted < len(sublists):
            upper_bound_slice = min((sublists_submitted + slices_per_iteration), len(sublists))
            cur_slice_start_indices = start_indices[sublists_submitted:upper_bound_slice]
            cur_slice_sublists = sublists[sublists_submitted:upper_bound_slice]

            # generate paths and initialize molecules (so that if they fail, this can be covered)
            tmp_output_dirs, tmp_input_mae_paths, \
            tmp_output_sdf_paths = self._generate_temporary_input_output_files(cur_slice_start_indices,
                                                                               cur_slice_sublists)

            # call "token guard" method (only executed, if block is specified in the configuration), which will wait
            # with the execution if not enough tokens are available at the moment
            self._apply_token_guard()

            # run in parallel; wait for all subjobs to finish before proceeding
            processes = []
            for chunk_index in range(len(tmp_output_dirs)):
                p = multiprocessing.Process(target=self._dock_subjob, args=(tmp_input_mae_paths[chunk_index],
                                                                            tmp_output_sdf_paths[chunk_index],
                                                                            tmp_output_dirs[chunk_index],
                                                                            number_ligands_per_sublist))
                processes.append(p)
                p.start()
            for p in processes:
                p.join()

            # add the number of input sublists rather than the output temporary folders to account for cases where
            # entire sublists failed to produce an input structure
            sublists_submitted += len(cur_slice_sublists)

            # parse the resulting sdf files
            for path_sdf_results in tmp_output_sdf_paths:
                # this is a protection against the case where empty (file size == 0 bytes) files are generated due to
                # a failure during docking
                if not os.path.isfile(path_sdf_results) or os.path.getsize(path_sdf_results) == 0:
                    continue

                for molecule in Chem.SDMolSupplier(path_sdf_results, removeHs=False):
                    if molecule is None:
                        continue
                    cur_conformer_name = str(molecule.GetProp("_Name"))

                    # add molecule to the appropriate ligand
                    for ligand in self.ligands:
                        if ligand.get_identifier() == cur_conformer_name:
                            ligand.add_conformer(molecule)
                            break

            # clean-up
            for path in tmp_output_dirs:
                shutil.rmtree(path)
            self._log_docking_progress(number_done=sublists_submitted, number_total=number_sublists)

        # sort the conformers (best to worst) and update their names to contain the conformer id
        # -> <ligand_number>:<enumeration>:<conformer_number>
        for ligand in self.ligands:
            ligand.set_conformers(sorted(ligand.get_conformers(),
                                         key=lambda x: float(x.GetProp(_ROE.GLIDE_DOCKING_SCORE)), reverse=False))
            ligand.add_tags_to_conformers()

        # log any docking fails
        self._docking_fail_check()

        # generate docking results as dataframe
        result_parser = GlideResultParser(ligands=[ligand.get_clone() for ligand in self.ligands])
        self._df_results = result_parser.as_dataframe()

        # set docking flag
        self._docking_performed = True

    def _dock_subjob(self, mae_ligand_path, path_sdf_results, tmp_output_dir, chunk_size):
        keywords = self._all_keywords()

        # 1) add "LIGANDFILE" keyword to list of keywords: full path to "mae" formatted ligands
        keywords[_EE.GLIDE_LIGANDFILE] = mae_ligand_path

        # 2) write the keyword-input file for the "Glide" backend; write-out to temporary file
        glide_keywords_path = gen_temp_file(suffix=".in", dir=tmp_output_dir)
        _ = self._write_keywords_to_file(keywords=keywords, path=glide_keywords_path)

        # 3) wait / sleep until job is completed
        path_tmp_results = os.path.join(os.path.dirname(glide_keywords_path),
                                        "".join([os.path.splitext(os.path.basename(glide_keywords_path))[0],
                                                 _ROE.GLIDE_SDF_DEFAULT_EXTENSION]))
        path_tmp_log = os.path.join(os.path.dirname(glide_keywords_path),
                                    "".join([os.path.splitext(os.path.basename(glide_keywords_path))[0],
                                             _ROE.GLIDE_LOG]))

        # 4) execute the "Glide" backend; note, that the first argument is the path to the keyword input file
        #    note: if the number of cores has been set, overwrite "N_JOBS" and parallelize internally and also note
        #    that each subjob requires a license; instead start each with "N_JOBS" = 1
        arguments = [glide_keywords_path]
        flags = deepcopy(self.parameters.glide_flags)
        flags[_EE.GLIDE_NJOBS] = 1

        # -WAIT leads to issues at times: The process may not return properly (e.g. because of writing problems) and
        # then gets stuck. Workaround with waiting for file completion, so remove it if set.
        flags.pop(_EE.GLIDE_WAIT, None)

        for key in flags.keys():
            arguments.append(key)
            if flags[key] != "":
                arguments.append(flags[key])
        execution_result = self._Schrodinger_executor.execute(command=_EE.GLIDE,
                                                              arguments=arguments,
                                                              check=False,
                                                              location=os.path.dirname(glide_keywords_path))

        # 5) check return code (anything but '0' is bad) and add "stdout" to log file
        time_exceeded = False
        if execution_result.returncode != 0:
            msg = f"Could not dock with Glide, error message: {execution_result.stdout}."
            self._logger.log(msg, _LE.ERROR)
            self._print_log_file(path_tmp_log)
            raise DockingRunFailed(msg)
        else:
            if self._wait_until_file_generation(path=path_tmp_results,
                                                path_log=path_tmp_log,
                                                interval_sec=10,
                                                maximum_sec=self._get_time_limit_per_ligand()*chunk_size) is False:
                time_exceeded = True
                self._logger.log(f"Sublist docking for output file {path_tmp_results} exceeded time limit or failed, all these ligands are ignored in the final write-out. This could mean that none of them could be docked or a runtime error in Glide occured.",
                                 _LE.DEBUG)

        # 6) load the log-file (if generated) and check if all went well
        if any_in_file(path_tmp_log, _EE.GLIDE_LOG_SUCCESS_STRING) and time_exceeded is False:
            self._logger.log(f"Finished sublist (input: {mae_ligand_path}, output: {path_sdf_results}).", _LE.DEBUG)
        else:
            self._print_log_file(path_tmp_log)

        # 7) collect the results; Glide outputs the sdf with a given, semi-hard-coded path; extract the sdf file
        if os.path.isfile(path_tmp_results):
            with gzip.open(path_tmp_results, "rb") as fin:
                with open(path_sdf_results, "wb") as fout:
                    shutil.copyfileobj(fin, fout)

    def write_docked_ligands(self, path, mode="all"):
        """This method overrides the parent class, docker.py write_docked_ligands method. This method writes docked
        ligands binding poses and conformers to a file. There is the option to output the best predicted binding pose
        per ligand, the best predicted binding pose per enumeration, or all the predicted binding poses

        :param path: Contains information on results output path
        :type path: string
        :param mode: Determines whether the output contains the best predicted binding pose per ligand, the best
            predicted binding pose per enumeration, or all the predicted binding poses
        :type mode: string, optional, default value is "all". Other possible values are "best_per_ligand" and
            "best_per_enumeration"
        :raises DockingRunFailed Error: This error is raised if the docking run has not been performed
        :raises OpenEye (OE) Fatal Error: This error is raised if the output file was unable to be created. Issues may
            be due to problems with the ligand structure
        :raises ValueError: This error is raised if the ligands are neither RDkit nor OpenEye readable
        """
        self._write_docked_ligands(path, mode, mol_type=_LP.TYPE_RDKIT)

    def _print_log_file(self, path):
        if os.path.isfile(path):
            with open(path, 'r') as log_file:
                log_file_raw = log_file.readlines()
                self._logger.log(f"Printing log file {path}:\n", _LE.DEBUG)
                for line in log_file_raw:
                    self._logger_blank.log(line.rstrip("\n"), _LE.DEBUG)
                self._logger_blank.log("", _LE.DEBUG)
                self._logger.log("--- End file", _LE.DEBUG)

    def _wait_until_file_generation(self, path, path_log=None, interval_sec=1, maximum_sec=None) -> bool:
        counter = 0
        while not os.path.exists(path):
            # wait for an interval
            time.sleep(interval_sec)
            counter = counter + 1

            # if a Glide logfile path has been specified, check, whether critical messages indicating an abort are there
            # note, that we return "True" to indicate that the "file generation" has nevertheless been completed
            if path_log is not None:
                if any_in_file(path_log, _EE.GLIDE_LOG_FAIL_STRINGS):
                    self._logger.log(f"A critical error occurred in sublist execution.", _LE.WARNING)
                    self._print_log_file(path_log)
                    return True
                if any_in_file(path_log, _EE.GLIDE_LOG_FINISHED_STRINGS):
                    # log file indicates job is done; give a bit of leeway to ensure the writing is done
                    time.sleep(3)
                    break

            # if there's time left, proceed
            if maximum_sec is not None and counter * interval_sec >= maximum_sec:
                break
        if os.path.exists(path):
            return True
        else:
            return False

    def _get_time_limit_per_ligand(self):
        # for "SP" method, it can be expected to that about 90 s / ligand is required at most; use a bit extra
        if self.parameters.time_limit_per_compound is not None:
            custom_time_limit = int(self.parameters.time_limit_per_compound)
            self._logger.log(f"Using custom time limit of {custom_time_limit} seconds per compound for Glide docking.",
                             _LE.DEBUG)
            return custom_time_limit
        else:
            self._logger.log("No custom time limit set for Glide execution per compound - will use default of 120 seconds.",
                             _LE.DEBUG)
            return 120
