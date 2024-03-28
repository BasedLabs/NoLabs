from rdkit import Chem
from rdkit import Geometry

import openeye.oechem as oechem

from dockstream.utils.smiles import to_smiles, to_mol


# This code is taken from Caitlin Bannan, Application Scienist at OpenEye (gist.github.com/bannanc)
# No license specified

def RDkitMolToOpenEyeMol(molecule: Chem.Mol, bySMILES: bool):
    """Creates an OpenEye molecule object that is identical to the input RDkit molecule."""

    # if "bySMILES" is True, simply do: rdkit molecule -> SMILE -> openeye molecule
    if bySMILES:
        OpenEye_molecule = oechem.OEMol()
        oechem.OESmilesToMol(OpenEye_molecule, to_smiles(molecule))
        return OpenEye_molecule

    # openeye stores bond orders as integers regardless of aromaticity
    # in order to properly extract these, we need to have the "Kekulized" version of the rdkit mol
    kekul_mol = Chem.Mol(molecule)
    Chem.Kekulize(kekul_mol, True)

    oemol = oechem.OEMol()
    map_atoms = dict()  # {rd_idx: oe_atom}

    # setting chirality in openey requires using neighbor atoms
    # therefore we can't do it until after the atoms and bonds are all added
    chiral_atoms = dict()  # {rd_idx: openeye chirality}
    for rda in molecule.GetAtoms():
        rd_idx = rda.GetIdx()

        # create a new atom
        oe_a = oemol.NewAtom(rda.GetAtomicNum())
        map_atoms[rd_idx] = oe_a
        oe_a.SetFormalCharge(rda.GetFormalCharge())
        oe_a.SetAromatic(rda.GetIsAromatic())

        # If chiral, store the chirality to be set later
        tag = rda.GetChiralTag()
        if tag == Chem.CHI_TETRAHEDRAL_CCW:
            chiral_atoms[rd_idx] = oechem.OECIPAtomStereo_R
        if tag == Chem.CHI_TETRAHEDRAL_CW:
            chiral_atoms[rd_idx] = oechem.OECIPAtomStereo_S

    # Similar to chirality, stereochemistry of bonds in OE is set relative to their neighbors
    stereo_bonds = list()
    # stereo_bonds stores tuples in the form (oe_bond, rd_idx1, rd_idx2, OE stereo specification)
    # where rd_idx1 and 2 are the atoms on the outside of the bond
    # i.e. Cl and F in the example above
    aro_bond = 0
    for rdb in molecule.GetBonds():
        a1 = rdb.GetBeginAtomIdx()
        a2 = rdb.GetEndAtomIdx()

        # create a new bond
        newbond = oemol.NewBond(map_atoms[a1], map_atoms[a2])

        order = rdb.GetBondTypeAsDouble()
        if order == 1.5:
            # get the bond order for this bond in the kekulized molecule
            order = kekul_mol.GetBondWithIdx(rdb.GetIdx()).GetBondTypeAsDouble()
            newbond.SetAromatic(True)
        else:
            newbond.SetAromatic(False)
        newbond.SetOrder(int(order))

        # determine if stereochemistry is needed
        tag = rdb.GetStereo()
        if tag == Chem.BondStereo.STEREOCIS or tag == Chem.BondStereo.STEREOZ:
            stereo_atoms = rdb.GetStereoAtoms()
            stereo_bonds.append((newbond, stereo_atoms[0], stereo_atoms[1], oechem.OEBondStereo_Cis))

            bond2 = molecule.GetBondBetweenAtoms(stereo_atoms[0], a1)
            bond4 = molecule.GetBondBetweenAtoms(stereo_atoms[1], a2)
        if tag == Chem.BondStereo.STEREOTRANS or tag == Chem.BondStereo.STEREOE:
            stereo_atoms = rdb.GetStereoAtoms()
            stereo_bonds.append((newbond, stereo_atoms[0], stereo_atoms[1], oechem.OEBondStereo_Trans))
            bond2 = molecule.GetBondBetweenAtoms(stereo_atoms[0], a1)
            bond4 = molecule.GetBondBetweenAtoms(stereo_atoms[1], a2)

    # Now that all of the atoms are connected we can set stereochemistry
    # starting with atom chirality
    for rd_idx, chirality in chiral_atoms.items():
        # chirality is set relative to neighbors, so we will get neighboring atoms
        # assign Right handed direction, check the cip stereochemistry
        # if the cip stereochemistry isn't correct then we'll set left and double check

        oea = map_atoms[rd_idx]
        neighs = [n for n in oea.GetAtoms()]
        # incase you look at the documentation oe has two options for handedness for example:
        # oechem.OEAtomStereo_Left == oechem.OEAtomStereo_LeftHanded
        oea.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Right)
        cip = oechem.OEPerceiveCIPStereo(oemol, oea)
        if cip != chirality:
            oea.SetStereo(neighs, oechem.OEAtomStereo_Tetra, oechem.OEAtomStereo_Left)
            new_cip = oechem.OEPerceiveCIPStereo(oemol, oea)
            if new_cip != chirality:
                # Note, I haven't seen this happen yet, but it shouldn't be a problem since there
                # is only 2 directions for handedness and we're only running this for chiral atoms
                pass

    # Set stereochemistry using the reference atoms extracted above
    for oeb, idx1, idx2, oestereo in stereo_bonds:
        oeb.SetStereo([map_atoms[idx1], map_atoms[idx2]], oechem.OEBondStereo_CisTrans, oestereo)

    # If the rdmol has a conformer, add its coordinates to the oemol
    # Note, this currently only adds the first conformer, it will need to be adjusted if the
    # you wanted to convert multiple sets of coordinates
    if molecule.GetConformers():
        conf = molecule.GetConformer()
        for rd_idx, oeatom in map_atoms.items():
            coords = conf.GetAtomPosition(rd_idx)
            oemol.SetCoords(oeatom, oechem.OEFloatArray(coords))

    # If RDMol has a title save it
    if molecule.HasProp("_Name"):
        oemol.SetTitle(molecule.GetProp("_Name"))

    # Clean Up phase
    # The only feature of a molecule that wasn't perceived above seemed to be ring connectivity, better to run it
    # here then for someone to inquire about ring sizes and get 0 when it shouldn't be
    oechem.OEFindRingAtomsAndBonds(oemol)

    # Fix the dimensionality (got the coordinates but not the appropriate tag - this is introduces with 2.0.something)
    oechem.OESetDimensionFromCoords(oemol)

    return oemol


def OpenEyeMolToRDkitMol(molecule: oechem.OEMol, bySMILES: bool):
    """Creates an OpenEye molecule object that is identical to the input RDkit molecule."""

    # if "bySMILES" is True, simply do: openeye molecule -> SMILE -> rdkit molecule
    if bySMILES:
        return to_mol(oechem.OEMolToSmiles(molecule))

    rdmol = Chem.RWMol()

    # RDKit keeps bond order as a type instead using these values, I don't really understand 7,
    # I took them from Shuzhe's example linked above
    _bondtypes = {1: Chem.BondType.SINGLE,
                  1.5: Chem.BondType.AROMATIC,
                  2: Chem.BondType.DOUBLE,
                  3: Chem.BondType.TRIPLE,
                  4: Chem.BondType.QUADRUPLE,
                  5: Chem.BondType.QUINTUPLE,
                  6: Chem.BondType.HEXTUPLE,
                  7: Chem.BondType.ONEANDAHALF, }

    # atom map lets you find atoms again
    map_atoms = dict()  # {oe_idx: rd_idx}
    for oea in molecule.GetAtoms():
        oe_idx = oea.GetIdx()
        rda = Chem.Atom(oea.GetAtomicNum())
        rda.SetFormalCharge(oea.GetFormalCharge())
        rda.SetIsAromatic(oea.IsAromatic())

        # unlike OE, RDK lets you set chirality directly
        cip = oechem.OEPerceiveCIPStereo(molecule, oea)
        if cip == oechem.OECIPAtomStereo_S:
            rda.SetChiralTag(Chem.CHI_TETRAHEDRAL_CW)
        if cip == oechem.OECIPAtomStereo_R:
            rda.SetChiralTag(Chem.CHI_TETRAHEDRAL_CCW)

        map_atoms[oe_idx] = rdmol.AddAtom(rda)

    # As discussed above, setting bond stereochemistry requires neighboring bonds
    # so we will store that information by atom index in this list
    stereo_bonds = list()
    # stereo_bonds will have tuples with the form (rda1, rda2, rda3, rda4, is_cis)
    # where rda[n] is an atom index for a double bond of form 1-2=3-4
    # and is_cis is a Boolean is True then onds 1-2 and 3-4 are cis to each other

    aro_bond = 0
    for oeb in molecule.GetBonds():
        # get neighboring rd atoms
        rd_a1 = map_atoms[oeb.GetBgnIdx()]
        rd_a2 = map_atoms[oeb.GetEndIdx()]

        # AddBond returns the total number of bonds, so addbond and then get it
        rdmol.AddBond(rd_a1, rd_a2)
        rdbond = rdmol.GetBondBetweenAtoms(rd_a1, rd_a2)

        # Assign bond type, which is based on order unless it is aromatic
        order = oeb.GetOrder()
        if oeb.IsAromatic():
            rdbond.SetBondType(_bondtypes[1.5])
            rdbond.SetIsAromatic(True)
        else:
            rdbond.SetBondType(_bondtypes[order])
            rdbond.SetIsAromatic(False)

        # If the bond has specified stereo add the required information to stereo_bonds
        if oeb.HasStereoSpecified(oechem.OEBondStereo_CisTrans):
            # OpenEye determined stereo based on neighboring atoms so get two outside atoms
            n1 = [n for n in oeb.GetBgn().GetAtoms() if n != oeb.GetEnd()][0]
            n2 = [n for n in oeb.GetEnd().GetAtoms() if n != oeb.GetBgn()][0]

            rd_n1 = map_atoms[n1.GetIdx()]
            rd_n2 = map_atoms[n2.GetIdx()]

            stereo = oeb.GetStereo([n1, n2], oechem.OEBondStereo_CisTrans)
            if stereo == oechem.OEBondStereo_Cis:
                stereo_bonds.append((rd_n1, rd_a1, rd_a2, rd_n2, True))
            elif stereo == oechem.OEBondStereo_Trans:
                stereo_bonds.append((rd_n1, rd_a1, rd_a2, rd_n2, False))

                # add bond stereochemistry:
    for (rda1, rda2, rda3, rda4, is_cis) in stereo_bonds:
        # get neighbor bonds
        bond1 = rdmol.GetBondBetweenAtoms(rda1, rda2)
        bond2 = rdmol.GetBondBetweenAtoms(rda3, rda4)

        # Since this is relative, the first bond always goes up
        # as explained above these names come from SMILES slashes so UP/UP is Trans and Up/Down is cis
        bond1.SetBondDir(Chem.BondDir.ENDUPRIGHT)
        if is_cis:
            bond2.SetBondDir(Chem.BondDir.ENDDOWNRIGHT)
        else:
            bond2.SetBondDir(Chem.BondDir.ENDUPRIGHT)

    # if oemol has coordinates (The dimension is non-zero)
    # add those coordinates to the rdmol
    if oechem.OEGetDimensionFromCoords(molecule) > 0:
        conformer = Chem.Conformer()
        oecoords = molecule.GetCoords()
        for oe_idx, rd_idx in map_atoms.items():
            (x, y, z) = oecoords[oe_idx]
            conformer.SetAtomPosition(rd_idx, Geometry.Point3D(x, y, z))
        rdmol.AddConformer(conformer)

        # Save the molecule title
    rdmol.SetProp("_Name", molecule.GetTitle())

    # Cleanup the rdmol
    # Note I copied UpdatePropertyCache and GetSSSR from Shuzhe's code to convert oemol to rdmol here:
    rdmol.UpdatePropertyCache(strict=False)
    Chem.GetSSSR(rdmol)
    # I added AssignStereochemistry which takes the directions of the bond set
    # and assigns the stereochemistry tags on the double bonds
    Chem.AssignStereochemistry(rdmol, force=False)

    # translate tags
    for pair in oechem.OEGetSDDataPairs(molecule):
        rdmol.SetProp(pair.GetTag(), pair.GetValue())

    return rdmol.GetMol()
