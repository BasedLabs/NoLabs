import multiprocessing
from typing import Optional, Any, List

import openeye.oeomega as oeomega
import openeye.oechem as oechem
import openeye.oedocking as oedocking
from copy import deepcopy

from pydantic import BaseModel
from typing_extensions import Literal

from dockstream.core.Schrodinger.Glide_docker import Parallelization
from dockstream.core.docker import Docker
from dockstream.core.OpenEye.OpenEye_result_parser import OpenEyeResultParser
from dockstream.utils.enums.OpenEye_enums import OpenEyeDockingConfigurationEnum, OpenEyeLigandPreparationEnum
from dockstream.utils.enums.docking_enum import DockingConfigurationEnum
from dockstream.utils.enums.logging_enums import LoggingConfigEnum

from dockstream.utils.translations.molecule_translator import MoleculeTranslator
from dockstream.utils.parallelization.general_utils import split_into_sublists
from dockstream.utils.dockstream_exceptions import DockingRunFailed

_LE = LoggingConfigEnum()
_LP = OpenEyeLigandPreparationEnum()
_CE = OpenEyeDockingConfigurationEnum()
_DE = DockingConfigurationEnum()


class OpenEyeParameters(BaseModel):
    prefix_execution: Optional[str] = None
    binary_location: Optional[str] = None
    parallelization: Optional[Parallelization]
    receptor_paths: List[str]
    scoring: str
    resolution: str
    number_poses: int

    def get(self, key: str) -> Any:
        """Temporary method to support nested_get"""
        return self.dict()[key]


class OpenEye(Docker):
    """Interface to the OpenEye backend."""

    backend: Literal["OpenEye"] = "OpenEye"
    parameters: OpenEyeParameters

    _docker: oedocking.OEDock = None
    _omega: oeomega.OEOmega = None

    class Config:
        underscore_attrs_are_private = True

    def __init__(self, **run_parameters):
        super().__init__(**run_parameters)
        self._logger.log("Using the OpenEye API docking version is deprecated - use Hybird.", _LE.WARNING)

    def _initialize_docker(self):
        # set the options for the OpenEye backend
        omegaOpts = oeomega.OEOmegaOptions()
        omegaOpts.SetStrictStereo(False)
        self._omega = oeomega.OEOmega(omegaOpts)
        oechem.OEThrow.SetLevel(10000)

        # load the target
        oereceptor = oechem.OEGraphMol()
        oedocking.OEReadReceptorFile(oereceptor, self.parameters.receptor_paths[0])

        # initialize the docker
        scoring_function_value = self._get_scoring_function_value(method=self.parameters.scoring)
        resolution_value = self._get_resolution_value(resolution=self.parameters.resolution)
        self._docker = oedocking.OEDock(scoring_function_value, resolution_value)
        self._docker.Initialize(oereceptor)
        self._logger.log(f"Initialized OpenEye receptor and set parameters.", _LE.DEBUG)

    def _get_score_from_conformer(self, conformer):
        return float(conformer.GetEnergy())

    def _get_scoring_function_value(self, method: str) -> int:
        # get the integer value that corresponds to the scoring function (e.g. "Hybrid1" == 18) for the initialization
        # notes: (1) there is also a "OEHybrid" version of "OEDock", but this differs only in the default setting for
        #            the scoring function
        #        (2) any string not matching one of the defined scoring functions, returns the same large integer, which
        #            is used for checks
        value = oedocking.OEDockMethodGetValue(method)
        if value == _CE.SCORING_INVALID_VALUE:
            raise DockingRunFailed("Scoring function " + method + " not supported.")
        return value

    def _get_resolution_value(self, resolution: str) -> int:
        value = oedocking.OESearchResolutionGetValue(resolution)
        if value == _CE.RESOLUTION_INVALID_VALUE:
            raise DockingRunFailed("Specified resolution " + resolution + " not supported.")
        return value

    def add_molecules(self, molecules: list):
        """This method overrides the parent class, docker.py add_molecules method. This method appends prepared
        ligands to a list for subsequent docking

        :param molecules: A list that is to contain all prepared ligands for subsequent docking
        :type molecules: list
        :raises NotImplementedError: Each backend must override the parent class, docker.py add_molecules method.
            Inability to do so or a bug causing incorrect implementation will raise a NotImplementedError
        """
        mol_trans = MoleculeTranslator(self.ligands, force_mol_type=_LP.TYPE_OPENEYE)
        mol_trans.add_molecules(molecules)
        self.ligands = mol_trans.get_as_openeye()
        self._docking_performed = False

    def get_sublists_for_docking(self, number_cores):
        """This method overrides the paren class, docker.py get_sublists_for_docking method. This method splits the
        ligands into sublists for docking to take advantage of parallel computing using >1 processing core on your
        computer. An additional caveat is that OpenEye cannot handle more than 10 ligands per sublist

        :param number_cores: Number of cores you want to allocate to the docking job. Ensure at least 1 core
            is kept free so other tasks can still run on your computer (ex. if you have 8 cores, allocate at most 7)
        :type number_cores: int
        :return: split_into_sublists containing ligands split into sublists for subsequent parallel docking
        """
        # TODO: fix the way OpenEye is parallelized to make handling consistent
        # decide how to slice the ligand list depending on whether a maximum length is defined or not
        slice_size = min(10, len(self.ligands))
        if self.parameters.parallelization is not None and \
            self.parameters.parallelization.max_compounds_per_subjob is not None:
            slice_size = min(max(self.parameters.parallelization.max_compounds_per_subjob, 1),
                             max(slice_size, 1))
        return split_into_sublists(input_list=self.ligands, partitions=None, slice_size=slice_size)

    def _dock(self, number_cores):
        # this function is quite complicated: conformers in OpenEye cannot be serialized (pickled), so we need to
        # wrap them in molecules instead and retrieve them individually afterwards

        # note, that the context for each child process is duplicated, so use a "Queue" to collect the results
        # submit the subjobs in batches of a certain maximum size, as "Queue" cannot handle bigger stuff; slice_size of
        # ten should be on the safe side for most applications, maybe review that part in the future or make parameter

        self._initialize_docker()

        _, sublists = self.get_sublists_for_docking(number_cores)
        number_sublists = len(sublists)
        self._logger.log(f"Split ligands into {len(sublists)} sublists for docking.", _LE.DEBUG)

        sublists_submitted = 0
        while sublists_submitted < len(sublists):
            processes = []
            return_queues = []
            for _ in range(number_cores):
                if sublists_submitted >= len(sublists):
                    continue
                cur_queue = multiprocessing.Queue()
                p = multiprocessing.Process(target=self._dock_subjob, args=(sublists[sublists_submitted],
                                                                            cur_queue))
                processes.append(p)
                p.start()
                return_queues.append(cur_queue)
                sublists_submitted += 1
            for p in processes:
                p.join()

            for cur_slice in [q.get() for q in return_queues]:
                for cur_ligand_name in cur_slice.keys():
                    for ligand in self.ligands:
                        if cur_ligand_name == ligand.get_identifier():
                            ligand.set_conformers(cur_slice[cur_ligand_name])
                            break
            self._log_docking_progress(number_done=sublists_submitted, number_total=number_sublists)

        # update conformer names to contain the conformer id -> <ligand_number>:<enumeration>:<conformer_number>
        for ligand in self.ligands:
            ligand.add_tags_to_conformers()

        # use self.docker.GetHighScoresAreBetter() to get a boolean indicating whether the particular scoring function
        # used considers "lower" values to be better or not; but the compounds are ordered properly here
        result_parser = OpenEyeResultParser(ligands=[ligand.get_clone() for ligand in self.ligands])
        self._df_results = result_parser.as_dataframe()
        self._logger.log(f"Successfully docked {len(self.ligands)} molecules.", _LE.DEBUG)

        # set the docking flag
        self._docking_performed = True

    def _dock_subjob(self, sublist, return_queue):
        dict_conformers = {}
        for ligand in sublist:
            mol = deepcopy(ligand.get_molecule())
            if mol is not None and self._omega(mol):
                dict_conformers[ligand.get_identifier()] = []
                docked_mol = oechem.OEMol()
                return_code = self._docker.DockMultiConformerMolecule(docked_mol,
                                                                      mol,
                                                                      self.parameters.number_poses)

                # add the score for each pose to the tag name of the scoring function, e.g. "Chemgauss3" and add each
                # pose for the ligands
                for conf in docked_mol.GetConfs():
                    # as only full-fledged "OEMol"s are serializable, make it a full-fledged molecule (a "GraphMol", as
                    # it contains only one conformation)
                    conf = oechem.OEGraphMol(conf)
                    oedocking.OESetSDScore(conf, self._docker, self.parameters.scoring)
                    dict_conformers[ligand.get_identifier()].append(conf)

        return_queue.put(dict_conformers)
        self._logger.log(f"Finished sublist with {len(sublist)} input molecules.", _LE.DEBUG)

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
        self._write_docked_ligands(path, mode, mol_type=_LP.TYPE_OPENEYE)

    def _load_ligands(self, path) -> list:
        ligands = []
        ifs = oechem.oemolistream()
        if not ifs.open(path):
            oechem.OEThrow.Fatal("Unable to open file for reading.")
        for mol in ifs.GetOEMols():
            ligands.append(mol)
        return ligands
