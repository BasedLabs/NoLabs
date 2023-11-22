import glob
import itertools
import subprocess
import os

from openmm.app import *
from openmm import *
from openmm.unit import *
from werkzeug.datastructures import FileStorage

from src.server.initializers.loggers import logger

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

pdb_tmp_file = 'tmp.pdb'
output_tmp_file = 'output.pdb'


def remove_conformations_backups():
    files = glob.glob('*.rtp*')
    for f in files:
        os.remove(f)
    files = glob.glob('*.itp*')
    for f in files:
        os.remove(f)
    files = glob.glob('*.gro*')
    for f in files:
        os.remove(f)


def remove_pdbs():
    files = glob.glob('*.pdb*')
    for f in files:
        os.remove(f)


amber_files = {
    'protein': [
        'amber14-all.xml',
        'amber14/protein.ff14SB.xml',
        'amber14/protein.ff15ipq.xml',
    ],
    'dna': [
        'amber14/DNA.OL15.xml',
        'amber14/DNA.bsc1.xml',
    ],
    'rna': [
        'amber14/RNA.OL3.xml',
    ],
    'lipid': [
        'amber14/lipid17.xml',
    ],
    'carbohydrates': [
        'amber14/GLYCAM_06j-1.xml',
    ],
    'water': [
        'amber14/tip3p.xml',
        'amber14/tip3pfb.xml',
        'amber14/tip4pew.xml',
        'amber14/tip4pfb.xml',
        'amber14/spce.xml',
        'amber14/opc.xml',
        'amber14/opc3.xml'
    ]
}

charmm = {
    'protein': ['charmm36.xml'],
    'water': [
        'charmm36/water.xml',
        'charmm36/spce.xml',
        'charmm36/tip3p-pme-b.xml',
        'charmm36/tip3p-pme-f.xml',
        'charmm36/tip4pew.xml',
        'charmm36/tip4p2005.xml',
        'charmm36/tip5p.xml',
        'charmm36/tip5pew.xml'
    ]
}

gromacs_force_fields = {
    'protein': ['amber03', 'amber94', 'amber96', 'amber99', 'amber99sb-ildn', 'amber99sb', 'amberGS', 'charmm27'],
    'water': ['spc', 'spce', 'ip3p', 'tip4p', 'tip5p', 'tips3p']
}


def run_simulation_on_pdb(input_pdb, output_pdb, force_fields, report_every=1000, report_steps=10000):
    try:
        pdb = PDBFile(input_pdb)
        forcefield = ForceField(*force_fields)
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1 * nanometer,
                                         constraints=HBonds)
        integrator = LangevinIntegrator(1 * kelvin, 1 / picosecond, 0.002 * picoseconds)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter(output_pdb, report_every))
        simulation.step(report_steps)
        logger.info("Success!")
        return True
    except Exception as e:
        logger.error(e)
        return False


def run_simulation_on_gromacs(input_top, input_gro, output_pdb, report_every=1000, report_steps=10000):
    try:
        gro = GromacsGroFile(input_gro)
        top = GromacsTopFile(input_top, periodicBoxVectors=gro.getPeriodicBoxVectors(),
                             includeDir='/usr/local/gromacs/share/gromacs/top')
        system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * nanometer,
                                  constraints=HBonds)
        integrator = LangevinMiddleIntegrator(1 * kelvin, 1 / picosecond, 0.002 * picoseconds)
        simulation = Simulation(top.topology, system, integrator)
        simulation.context.setPositions(gro.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter(output_pdb, report_every))
        simulation.step(report_steps)
        return True
    except Exception as e:
        logger.error(e)
        return False


def fix_pdb(input_pdb):
    import pdbfixer
    import openmm.app as mm_app

    output_pdb = f'{input_pdb}-fixed.pdb'
    if os.path.exists(output_pdb):
        return output_pdb

    fixer = pdbfixer.PDBFixer(input_pdb)
    logger.info("Finding missing residues")
    #fixer.findMissingResidues()
    logger.info("Finding nonstandard residues")
    #fixer.findNonstandardResidues()
    logger.info('Replacing nonstandard residues')
    #fixer.replaceNonstandardResidues()
    logger.info('Finding missing atoms')
    fixer.findMissingAtoms()
    logger.info('Adding missing atoms')
    fixer.addMissingAtoms()
    logger.info('Adding missing hydrogens')
    fixer.addMissingHydrogens(7.0)
    logger.info('Writing fixed file')
    mm_app.PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))
    return output_pdb


def generate_gromacs_files(input_pdb, output_gro, topology_top, force_field, water):
    res = subprocess.run(['/usr/local/gromacs/bin/gmx', 'pdb2gmx', '-f',
                          input_pdb, '-o',
                          output_gro, '-p',
                          topology_top,
                          '-ignh', '-ff', force_field.lower(), '-water', water.lower()],
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
    if os.path.exists(output_gro) and os.path.exists(topology_top):
        logger.info(res.stdout.decode('utf-8'))
        return True
    logger.error(res.stderr.decode('utf-8'))
    return False


def permute_simulation(pdb_content: FileStorage, report_every=1000, simulations_count=10000):
    remove_conformations_backups()
    remove_pdbs()
    pdb_content.save(pdb_tmp_file)

    def inner():
        tries = [amber_files, charmm]

        for input_pdb in [fix_pdb(pdb_tmp_file), pdb_tmp_file]:
            #for water in gromacs_force_fields['water']:
            #    for force_field in gromacs_force_fields['protein']:
            #        logger.info(
            #            f'Trying to process {input_pdb} by GROMACS using force field {force_field} and water {water}')
            #        gro = f'{input_pdb}.gro'
            #        top = f'{input_pdb}.top'
            #        res = generate_gromacs_files(input_pdb, f'{input_pdb}.gro', f'{input_pdb}.top', force_field,
            #                          water)
            #        if res and run_simulation_on_gromacs(top, gro, output_tmp_file, report_every, simulations_count):
            #            logger.info('Success')
            #            return
            #        remove_conformations_backups()

            for t in tries:
                for protein_force_fields in [pff for pff in t['protein']]:
                    for water_force_fields in [wff for wff in t['water']]:
                        logger.info(
                            f'Trying to process {input_pdb} by OPENMM using force field {protein_force_fields} and water {water_force_fields}')

                        sim_res = run_simulation_on_pdb(input_pdb, output_tmp_file,
                                                        [protein_force_fields, water_force_fields])
                        if sim_res:
                            logger.info('Success')
                            return
                        remove_conformations_backups()

                        other_force_fields_keys = [ff_key for ff_key in
                                                   [key for key in t.keys() if key != 'water' and key != 'protein']]

                        for key in other_force_fields_keys:
                            force_fields = itertools.chain([protein_force_fields, water_force_fields], t[key])
                            sim_res = run_simulation_on_pdb(input_pdb,
                                                            output_tmp_file,
                                                            force_fields, report_every, simulations_count)
                            if sim_res:
                                logger.info('Success')
                                return
                            remove_conformations_backups()
    try:
        inner()

        if os.path.exists(output_tmp_file):
            with open(output_tmp_file, 'r') as f:
                pdb_output = f.read()
                return pdb_output
        else:
            logger.info('Unable to generate conformations for this .pdb file')
            return
    finally:
        remove_conformations_backups()
        remove_pdbs()
