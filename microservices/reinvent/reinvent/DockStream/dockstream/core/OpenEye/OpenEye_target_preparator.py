import os

from dockstream.core.target_preparator import TargetPreparator
import openeye.oechem as oechem
import openeye.oedocking as oedocking

from dockstream.utils.dockstream_exceptions import TargetPreparationFailed

from dockstream.utils.enums.OpenEye_enums import OpenEyeTargetPreparationEnum
from dockstream.containers.target_preparation_container import TargetPreparationContainer


class OpenEyeTargetPreparator(TargetPreparator):
    """Class that deals with all the target preparatory steps needed before docking using "OpenEye" can commence."""

    def __init__(self, conf: TargetPreparationContainer, target, run_number=0):
        self._TP = OpenEyeTargetPreparationEnum()

        # invoke base class's constructor first
        super().__init__(conf=conf, run_number=run_number)

        # check, whether the backend run specified is an "OpenEye" one
        if self._run_parameters[self._TP.RUNS_BACKEND] != self._TP.RUNS_BACKEND_OPENEYE:
            raise TargetPreparationFailed("Tried to make an OpenEye preparation with different backend specification.")

        # treat the target: either load a file or store the molecule internally
        if isinstance(target, str):
            if os.path.isfile(target):
                _, file_extension = os.path.splitext(target)
                if file_extension == ".pdb":
                    istream = oechem.oemolistream(target)
                    protein = oechem.OEGraphMol()
                    oechem.OEReadMolecule(istream, protein)
                    self._protein = protein
                else:
                    raise TargetPreparationFailed("Specified input file must be in PDB format for OpenEye.")
            else:
                raise TargetPreparationFailed("Input target file does not exist.")
        elif isinstance(target, oechem.OEGraphMol):
            self._target = target
        else:
            raise TargetPreparationFailed("Constructor only accepts an OEGraphMol (OpenEye) object or a file path.")
        self._logger.log("Stored target as OpenEye molecule.", self._TL.DEBUG)

    def specify_cavity(self):
        target = oechem.OEGraphMol()
        if self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD] == self._TP.CAVITY_METHOD_BOX:
            # specify 6 floating-point numbers which define a "box" around the cavity
            limits = [float(x) for x in self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_BOX_LIMITS]]
            if len(limits) != 6:
                raise TargetPreparationFailed("The limits for the box specification must be an array of 6 values.")
            box = oedocking.OEBox(*limits)
            oedocking.OEMakeReceptor(target, self._protein, box)
        elif self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD] == self._TP.CAVITY_METHOD_REFERENCE:
            # use a reference molecule to specify the cavity
            ref_istream = oechem.oemolistream(self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_PATH])

            # set the provided format
            ref_format = self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_REFERENCE_FORMAT].upper()
            if ref_format == self._TP.CAVITY_REFERENCE_FORMAT_SDF:
                ref_istream.SetFormat(oechem.OEFormat_SDF)
            elif ref_format == self._TP.CAVITY_REFERENCE_FORMAT_PDB:
                ref_istream.SetFormat(oechem.OEFormat_PDB)
            else:
                raise TargetPreparationFailed("Specified format not supported!")

            ref = oechem.OEGraphMol()
            oechem.OEReadMolecule(ref_istream, ref)
            oedocking.OEMakeReceptor(target, self._protein, ref)
        elif self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD] == self._TP.CAVITY_METHOD_HINT:
            # specify a point (a "hint") that is in or close at the cavity; three coordinates required
            coordinates = [float(x) for x in self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_HINT_COORDINATES]]
            if len(coordinates) != 3:
                raise TargetPreparationFailed("Method hint requires 3 values.")
            oedocking.OEMakeReceptor(target, self._protein, *coordinates)
        else:
            raise TargetPreparationFailed("Specified cavity determination method not defined for OpenEye.")
        self._target = target
        self._logger.log(f"Generated OpenEye receptor with method {self._run_parameters[self._TP.CAVITY][self._TP.CAVITY_METHOD]}.", self._TL.DEBUG)

    def write_target(self, path):
        _, file_extension = os.path.splitext(path)
        if file_extension != ".oeb" and file_extension != ".oeb.gz":
            raise TargetPreparationFailed("Receptor files must end on either .oeb or .oeb.gz.")
        oedocking.OEWriteReceptorFile(self._target, path)
        self._logger.log(f"Wrote receptor to file {path}.", self._TL.DEBUG)
