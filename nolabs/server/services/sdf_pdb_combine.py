import os.path

from rdkit import Chem
from rdkit.Chem import AllChem
import io

def combine_sdf_pdb(sdf_file_content, pdb_file_content):
    tmp_sdf = 'tmp_sdf.sdf'
    tmp_pdb = 'tmp_pdb.pdb'
    temp_res_pdb = 'temp_molecule.pdb'

    def cleanup():
        if os.path.exists(tmp_sdf):
            os.remove(tmp_sdf)

        if os.path.exists(tmp_pdb):
            os.remove(tmp_pdb)

        if os.path.exists(temp_res_pdb):
            os.remove(temp_res_pdb)

    cleanup()

    try:
        with open(tmp_sdf, 'w') as f:
            f.write(sdf_file_content)

        with open(tmp_pdb, 'w') as f:
            f.write(pdb_file_content)

        # Load SDF file
        suppl = Chem.SDMolSupplier(tmp_sdf)
        mol = next(suppl)  # Assuming only one molecule in the SDF
        if mol is None:
            raise ValueError("No molecule found in SDF file")

        # Generate 3D coordinates if not present
        if not mol.GetNumConformers():
            AllChem.EmbedMolecule(mol)
            AllChem.UFFOptimizeMolecule(mol)

        # Write the molecule from SDF to a temporary PDB file
        writer = Chem.PDBWriter(temp_res_pdb)
        writer.write(mol)
        writer.close()

        combined_file_content = ''
        # Combine the temporary PDB with the original PDB
        with open(temp_res_pdb, 'r') as temp_file, open(tmp_pdb, 'r') as original_file, io.StringIO() as combined_file:
            # Write the contents of the SDF-derived PDB
            for line in temp_file:
                combined_file.write(line)

            # Write the contents of the original PDB
            for line in original_file:
                combined_file.write(line)

            combined_file_content = combined_file.getvalue()
            return combined_file_content
    finally:
        cleanup()