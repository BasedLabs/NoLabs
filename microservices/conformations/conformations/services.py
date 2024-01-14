import os
import subprocess
import uuid

from conformations.integrator_factory import create_integrator
from conformations.api_models import RunPdbSimulationsRequest, RunSimulationsResponse, \
    RunGromacsSimulationsRequest, \
    RunPdbFixerRequest, RunPdbFixerResponse
from conformations.loggers import logger

from openmm.app import *
from openmm import *
import pdbfixer
import openmm.app as mm_app
from openmm.unit import *

from conformations.api_models import GenGroTopRequest, GenGroTopResponse

# TODO Add logging!


__all__ = ['run_pdb_fixer', 'run_pdb_simulation', 'run_gromacs_simulation', 'generate_gromacs_files']



def run_pdb_fixer(parameters: RunPdbFixerRequest) -> RunPdbFixerResponse:
    input_pdb = f'{uuid.uuid4()}.pdb'
    output_pdb = f'{uuid.uuid4()}.pdb'

    try:
        with open(input_pdb, 'w') as f:
            f.write(parameters.pdb_content)
        fixer = pdbfixer.PDBFixer(input_pdb)
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.findMissingAtoms()
        if parameters.replace_nonstandard_residues:
            fixer.replaceNonstandardResidues()
        if parameters.add_missing_atoms:
            fixer.addMissingAtoms()
        if parameters.add_missing_hydrogens:
            fixer.addMissingHydrogens(7.0)
        mm_app.PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))
        with open(output_pdb) as f:
            output_pdb_content = f.read()
        return RunPdbFixerResponse(pdb_content=output_pdb_content, errors=[])
    except Exception as e:
        logger.exception('Exception in pdb fixer')
        return RunPdbFixerResponse(pdb_content=None, errors=['Unable to run pdb fixer. Internal error', str(e)])
    finally:
        _remove_file_if_exists(output_pdb)
        _remove_file_if_exists(input_pdb)


def run_gromacs_simulation(parameters: RunGromacsSimulationsRequest) -> RunSimulationsResponse:
    output_pdb = f'{uuid.uuid4()}.pdb'
    input_gro = f'{uuid.uuid4()}.gro'
    input_top = f'{uuid.uuid4()}.top'

    friction_coeff = parameters.friction_coeff / picosecond
    step_size = parameters.step_size * picoseconds
    temperature_kelvin = parameters.temperature_k
    integrator_class = create_integrator(parameters.integrator)
    total_frames = parameters.total_frames
    take_frame_every = parameters.take_frame_every

    try:
        with open(input_gro, 'w') as f:
            f.write(parameters.gro)

        with open(input_top, 'w') as f:
            f.write(parameters.top)

        gro = GromacsGroFile(input_gro)
        top = GromacsTopFile(input_top, periodicBoxVectors=gro.getPeriodicBoxVectors(),
                             includeDir='/usr/local/gromacs/share/gromacs/top')
        system = top.createSystem(nonbondedMethod=PME, nonbondedCutoff=1 * nanometer,
                                  constraints=HBonds)
        integrator = integrator_class(temperature_kelvin, friction_coeff, step_size)
        simulation = Simulation(top.topology, system, integrator)
        simulation.context.setPositions(gro.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter(output_pdb, take_frame_every))
        simulation.step(total_frames)
        with open(output_pdb, 'r') as f:
            pdb_content = f.read()
        return RunSimulationsResponse(pdb_content=pdb_content, errors=[])
    except Exception as e:
        logger.exception('Exception in run gromacs simulations')
        return RunSimulationsResponse(pdb_content=None, errors=['Unable to generate simulations due to internal error', str(e)])
    finally:
        _remove_file_if_exists(output_pdb)
        _remove_file_if_exists(input_gro)
        _remove_file_if_exists(input_top)


def run_pdb_simulation(parameters: RunPdbSimulationsRequest) -> RunSimulationsResponse:
    output_pdb = f'{uuid.uuid4()}.pdb'
    input_pdb = f'{uuid.uuid4()}.pdb'

    friction_coeff = parameters.friction_coeff / picosecond
    step_size = parameters.step_size * picoseconds
    temperature_kelvin = parameters.temperature_k
    integrator_class = create_integrator(parameters.integrator)
    total_frames = parameters.total_frames
    take_frame_every = parameters.take_frame_every

    try:
        with open(input_pdb, 'w') as f:
            f.write(parameters.pdb_content)
        pdb = PDBFile(input_pdb)
        forcefield = ForceField(*[parameters.force_field.value, parameters.water_force_field.value])
        system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1 * nanometer,
                                         constraints=HBonds)
        integrator = integrator_class(temperature_kelvin, friction_coeff, step_size)
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter(output_pdb, take_frame_every))
        simulation.step(total_frames)
        with open(output_pdb, 'r') as f:
            pdb_content = f.read()
        return RunSimulationsResponse(pdb_content=pdb_content, errors=[])
    except Exception as e:
        logger.exception('Exception in pdb simulations')
        return RunSimulationsResponse(pdb_content=None, errors=['Unable to generate simulations due to internal error', str(e)])
    finally:
        _remove_file_if_exists(output_pdb)
        _remove_file_if_exists(input_pdb)


def generate_gromacs_files(parameters: GenGroTopRequest) -> GenGroTopResponse:
    input_pdb = f'{uuid.uuid4()}.pdb'
    output_gro = f'{uuid.uuid4()}.gro'
    output_top = f'{uuid.uuid4()}.top'

    try:
        with open(input_pdb, 'w') as f:
            f.write(parameters.pdb_content)
        program = ['/usr/local/gromacs/bin/gmx', 'pdb2gmx', '-f',
                   input_pdb, '-o',
                   output_gro, '-p',
                   output_top, '-ff', parameters.force_field.value, '-water', parameters.water_force_field.value]
        if parameters.ignore_missing_atoms:
            program.append('-ignh')
        res = subprocess.run(program, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        if os.path.exists(output_gro) and os.path.exists(output_top):
            with open(output_gro, 'r') as gro, open(output_top, 'r') as top:
                return GenGroTopResponse(gro=gro.read(), top=top.read())
        return GenGroTopResponse(gro=None, top=None, errors=[res.stderr.decode('utf-8')])
    except Exception as e:
        logger.exception('Exception in generate gromacs files')
        return GenGroTopResponse(gro=None, top=None, errors=['Unable to generate .gro and .top files from this pdb', str(e)])
    finally:
        _remove_file_if_exists(input_pdb)
        _remove_file_if_exists(output_gro)
        _remove_file_if_exists(output_top)


def _remove_file_if_exists(path: str) -> None:
    if os.path.exists(path):
        os.remove(path)