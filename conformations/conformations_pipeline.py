import glob
import itertools
import logging
import subprocess
from typing import Union
from urllib.parse import urljoin

import requests
from openmm.app import *
from openmm import *
from openmm.unit import *
import re

import settings

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

pdb_tmp_file = 'tmp.pdb'
output_tmp_file = 'output.pdb'

last_log = ''
last_error = ''


class LogsFilter(logging.Filter):
    def filter(self, record):
        global last_log
        if record.levelno == logging.INFO:
            logs_url = urljoin(settings.MAIN_SERVICE_URL, 'logging/conformations/logs')
            requests.post(logs_url, json={'get-logs': record.msg})
        return True


class ErrorsFilter(logging.Filter):
    def filter(self, record):
        global last_error
        if record.levelno == logging.ERROR:
            logs_url = urljoin(settings.MAIN_SERVICE_URL, 'logging/conformations/errors')
            requests.post(logs_url, json={'conformations-errors': record.msg})
        return True


logger = logging.getLogger('my_logger')
logger.setLevel(logging.DEBUG)

ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
ch.addFilter(LogsFilter())
logger.addHandler(ch)

conformations_errors_logger = logging.getLogger('conformations_logger')
conformations_errors_logger.setLevel(logging.DEBUG)

ch_conformations = logging.StreamHandler()
ch_conformations.setLevel(logging.DEBUG)
ch_conformations.addFilter(ErrorsFilter())
conformations_errors_logger.addHandler(ch_conformations)


class FileStorage:
    def __init__(self, data: str):
        self.data = data

    def save(self, path):
        with open(path, 'w') as f:
            f.write(self.data)


def pipeline(pdb_content: str,
             total_frames: int = 10000,
             take_frame_every: int = 1000,
             integrator: str = 'LangevinIntegator',
             friction_coeff: float = 1,
             step_size: float = 0.002,
             temperature_kelvin: float = 273.15,
             replace_nonstandard_residues: bool = True,
             add_missing_atoms: bool = True,
             add_missing_hydrogens: bool = True,
             ignore_missing: bool = False):
    '''
    Run simulation
    :param pdb_content: pdb content
    :param total_frames: total frames to take
    :param take_frame_every: take frame every nth
    :param integrator: integrator, one of ['LangevinIntegator', 'LangevinMiddleIntegrator', 'NoseHooverIntegrator', 'BrownianIntegrator', 'VariableVerletIntegrator'].
    See http://docs.openmm.org/7.1.0/api-python/library.html#integrators
    :param friction_coeff ( picoseconds, default is 1 / picosecond ): see http://docs.openmm.org/7.1.0/api-python/generated/simtk.openmm.openmm.LangevinIntegrator.html
    :param step_size: ( picoseconds, default is 0.002 / picoseconds ) see http://docs.openmm.org/7.1.0/api-python/generated/simtk.openmm.openmm.LangevinIntegrator.html
    :param replace_nonstandard_residues: whether to replace nonstandard residues that cannot be found in standard dbs
    :param add_missing_atoms: add missing atoms
    :param add_missing_hydrogens: add missing hydrogens
    :param ignore_missing: ignore missing. Dangerous. See https://manual.gromacs.org/documentation/current/onlinehelp/gmx-pdb2gmx.html -ignh option
    :return:
    '''
    global last_log
    global last_error

    last_log = ''
    last_error = ''

    pdb_content = FileStorage(pdb_content)
    meaningful_errors_cache = set()

    friction_coeff = friction_coeff / picosecond
    step_size = step_size * picoseconds

    def integrator_factory():
        if integrator == 'LangevinIntegator':
            return LangevinIntegrator
        if integrator == 'LangevinMiddleIntegrator':
            return LangevinMiddleIntegrator
        if integrator == 'NoseHooverIntegrator':
            return NoseHooverIntegrator
        if integrator == 'BrownianIntegrator':
            return BrownianIntegrator
        if integrator == 'VariableVerletIntegrator':
            return VariableVerletIntegrator

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

    def clean_logs():
        global last_log
        global last_error

        last_log = ''
        last_error = ''

    def log_error(e):
        if isinstance(e, Exception):
            error = str(e)
            meaningful_error = extract_meaningful_error(error)
            if meaningful_error:
                if meaningful_error not in meaningful_errors_cache:
                    meaningful_errors_cache.add(meaningful_error)
                    conformations_errors_logger.error(meaningful_error)
            else:
                logger.error(error[-100:])
        else:
            error = e.decode('utf-8')
            meaningful_error = extract_meaningful_error(error)
            if meaningful_error:
                if meaningful_error not in meaningful_errors_cache:
                    meaningful_errors_cache.add(meaningful_error)
                    conformations_errors_logger.error(meaningful_error)
            else:
                logger.error(error[:100])

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

    def extract_meaningful_error(error):
        r = re.search('No template.*?\.', error)
        if r:
            residue = r.group(0).replace('No template found for residue', '').strip().strip('.')
            return f"{residue} not found in force fields databases. Replace or fix it"
        r = re.search('.*?not found in force fields databases', error)  # 1 (HIS) not found in force fields databases
        if r:
            residue = r.group(0).replace('not found in force fields databases', '').strip().strip('.')
            return f"{residue} not found in force fields databases. Replace or fix it"

    def run_simulation_on_pdb(input_pdb, output_pdb, force_fields):
        try:
            pdb = PDBFile(input_pdb)
            forcefield = ForceField(*force_fields)
            system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1 * nanometer,
                                             constraints=HBonds)
            integrator_class = integrator_factory()
            integrator = integrator_class(temperature_kelvin, friction_coeff, step_size)
            simulation = Simulation(pdb.topology, system, integrator)
            simulation.context.setPositions(pdb.positions)
            simulation.minimizeEnergy()
            simulation.reporters.append(PDBReporter(output_pdb, take_frame_every))
            simulation.step(total_frames)
            logger.info("Success!")
            return True
        except Exception as e:
            log_error(e)
            return False

    def run_simulation_on_gromacs(input_top, input_gro, output_pdb):
        try:
            gro = GromacsGroFile(input_gro)
            top = GromacsTopFile(input_top, periodicBoxVectors=gro.getPeriodicBoxVectors(),
                                 includeDir='/usr/local/gromacs/share/gromacs/top')
            system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * nanometer,
                                      constraints=HBonds)
            integrator_class = integrator_factory()
            integrator = integrator_class(temperature_kelvin, friction_coeff, step_size)
            simulation = Simulation(top.topology, system, integrator)
            simulation.context.setPositions(gro.positions)
            simulation.minimizeEnergy()
            simulation.reporters.append(PDBReporter(output_pdb, take_frame_every))
            simulation.step(total_frames)
            return True
        except Exception as e:
            log_error(e)
            return False

    def fix_pdb(input_pdb):
        if replace_nonstandard_residues or add_missing_atoms or add_missing_hydrogens:
            import pdbfixer
            import openmm.app as mm_app
            output_pdb = f'{input_pdb}-fixed.pdb'
            if os.path.exists(output_pdb):
                return output_pdb
            fixer = pdbfixer.PDBFixer(input_pdb)
            logger.info("Finding missing residues")
            fixer.findMissingResidues()
            logger.info("Finding nonstandard residues")
            fixer.findNonstandardResidues()
            logger.info('Finding missing atoms')
            fixer.findMissingAtoms()
            if replace_nonstandard_residues:
                logger.info('Replacing nonstandard residues')
                fixer.replaceNonstandardResidues()
            if add_missing_atoms:
                logger.info('Adding missing atoms')
                fixer.addMissingAtoms()
            if add_missing_hydrogens:
                logger.info('Adding missing hydrogens')
                fixer.addMissingHydrogens(7.0)
            logger.info('Writing fixed file')
            mm_app.PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))
            return output_pdb
        return input_pdb

    def generate_gromacs_files(input_pdb, output_gro, topology_top, force_field, water):
        program = ['/usr/local/gromacs/bin/gmx', 'pdb2gmx', '-f',
                   input_pdb, '-o',
                   output_gro, '-p',
                   topology_top, '-ff', force_field.lower(), '-water', water.lower()]
        if ignore_missing:
            program.append('-ignh')
        res = subprocess.run(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        if os.path.exists(output_gro) and os.path.exists(topology_top):
            logger.info(res.stdout.decode('utf-8')[:100])
            return True
        log_error(res.stderr)
        return False

    def permute_simulation(pdb_content: FileStorage):
        remove_conformations_backups()
        remove_pdbs()
        pdb_content.save(pdb_tmp_file)

        def inner():
            tries = [amber_files, charmm]

            for input_pdb in [pdb_tmp_file, fix_pdb(pdb_tmp_file)]:
                for water in gromacs_force_fields['water']:
                    for force_field in gromacs_force_fields['protein']:
                        logger.info(
                            f'Trying to process your pdb by GROMACS using force field {force_field} and water {water}')
                        gro = f'{input_pdb}.gro'
                        top = f'{input_pdb}.top'
                        res = generate_gromacs_files(input_pdb, f'{input_pdb}.gro', f'{input_pdb}.top', force_field,
                                                     water)
                        if res and run_simulation_on_gromacs(top, gro, output_tmp_file):
                            logger.info('Success')
                            return
                        remove_conformations_backups()

                for t in tries:
                    for protein_force_fields in [pff for pff in t['protein']]:
                        for water_force_fields in [wff for wff in t['water']]:
                            logger.info(
                                f'Trying to process your pdb by OPENMM using force field {protein_force_fields} and water {water_force_fields}')

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
                                                                force_fields)
                                if sim_res:
                                    logger.info('Success')
                                    return
                                remove_conformations_backups()

        try:
            inner()

            if os.path.exists(output_tmp_file):
                with open(output_tmp_file, 'r') as f:
                    pdb_output = f.read()
                    conformations_errors_logger.info('<success>')
                    return pdb_output
            else:
                conformations_errors_logger.info('<end>')
                logger.error('Unable to generate conformations for this .pdb file')
                return
        finally:
            remove_conformations_backups()
            remove_pdbs()

    return permute_simulation(pdb_content)
