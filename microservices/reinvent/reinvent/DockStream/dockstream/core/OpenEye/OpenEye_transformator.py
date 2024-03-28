import openeye.oechem as oechem

from dockstream.utils.dockstream_exceptions import LigandPreparationFailed, TransformationFailed

from dockstream.utils.enums.OpenEye_enums import OpenEyeDockingConfigurationEnum, OpenEyeLigandPreparationEnum
from dockstream.core.transformator import Transformator


class OpenEyeTransformator(Transformator):
    """Class that applies SMIRKS (OpenEye style) to compounds before they are further processed."""

    def __init__(self, conf):
        super().__init__(conf)
        self._LP = OpenEyeLigandPreparationEnum()
        self._CE = OpenEyeDockingConfigurationEnum()

    def transform(self, ligands) -> list:
        number_input = len(ligands)
        failed_indices = []

        # code based on Graeme Robb's scrip
        if self._type == self._TE.TRANSFORMATION_TYPE_SMIRKS:
            rxn = oechem.OEUniMolecularRxn(self._smirk)

            for ligand_number, ligand in enumerate(ligands):
                try:
                    molecule = oechem.OEMol()
                    oechem.OESmilesToMol(molecule, ligand.get_smile())

                    # apply the smirk
                    oechem.OEAddExplicitHydrogens(molecule)
                    success = rxn(molecule)
                    if success:
                        ligand.set_smile(oechem.OEMolToSmiles(molecule))
                    else:
                        raise TransformationFailed
                except Exception as e:
                    failed_indices.append(ligand_number)

            if self._fail_action == self._TE.TRANSFORMATION_FAIL_ACTION_DISCARD:
                failed_indices.reverse()
                for index in failed_indices:
                    self._logger.log(f"Failed transformation, discarding ligand with smile {ligands[index].get_smile()}.", self._LE.DEBUG)
                    del ligands[index]
        else:
            self._logger.log(f"For transformation backend {self._backend}, only type {self._TE.TRANSFORMATION_TYPE_SMIRKS} is supported.", self._LE.ERROR)
            raise LigandPreparationFailed(f"For transformation backend {self._backend}, only type {self._TE.TRANSFORMATION_TYPE_SMIRKS} is supported.")

        self._logger.log(f"Of {number_input} input smiles, {len(ligands)} smiles were transformed / retained ({len(failed_indices)} transformations failed with \"fail_action\" set to {self._fail_action}).", self._LE.DEBUG)
        for ligand in ligands:
            self._logger_blank.log(ligand.get_smile(), self._LE.DEBUG)
        return ligands
