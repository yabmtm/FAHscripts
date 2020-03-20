#!/usr/bin/env python

from openforcefield.topology import Molecule
from openforcefield.typing.engines.smirnoff import ForceField
from simtk import openmm, unit
from simtk.openmm import app
from openmmforcefields.generators import SystemGenerator
import parmed
import pandas as pd
import sys, subprocess, os

water_model = 'tip3p'
solvent_padding = 10.0 * unit.angstrom
ionic_strength = 100 * unit.millimolar # 100
pressure = 1.0 * unit.atmospheres
collision_rate = 91.0 / unit.picoseconds
temperature = 300.0 * unit.kelvin
timestep = 2.0 * unit.femtoseconds
nsteps_equil = 50000 # test

protein_forcefield = 'amber14/protein.ff14SB.xml'

# trying to use SMIRNOFF causes a torsion error during SystemGenerator call
#small_molecule_forcefield = 'openff_unconstrained-1.1.0.offxml'
small_molecule_forcefield = 'openff-1.1.0'
#small_molecule_forcefield = 'gaff-2.11'
solvation_forcefield = 'amber14/tip3p.xml'

ligand_data = pd.read_pickle('../ligands.pkl')

for ligand_ndx in range(int(sys.argv[1]),int(sys.argv[2])): # range(1,41341):
  try:
    print(f'Processing RUN{ligand_ndx}')
    print(f'Receptor: {ligand_data.TransFSReceptor[ligand_ndx]}') 
    receptor_file = f'../fixed_{ligand_data.TransFSReceptor[ligand_ndx]}.pdb'

   # create the system, solvate, and write integrator
    output_prefix = f'Mpro-x0072_0_LIG{ligand_ndx}'
    if os.path.exists(f'{output_prefix}/conf.gro'):
        continue
    ligand_file = f'{output_prefix}/Mpro-x0072_0_LIG{ligand_ndx}_h.sdf'

    ligand = Molecule.from_file(ligand_file)

    receptor = app.PDBFile(receptor_file)
    receptor_structure = parmed.load_file(receptor_file)
    ligand_structure = parmed.load_file(f'{output_prefix}/Mpro-x0072_0_LIG{ligand_ndx}_h.pdb')
    complex_structure = receptor_structure + ligand_structure

    barostat = openmm.MonteCarloBarostat(pressure, temperature)
    forcefield_kwargs = {'removeCMMotion': False, 'ewaldErrorTolerance': 5e-04,
        'nonbondedMethod': app.PME, 'constraints': None, 'rigidWater': False}
    system_generator = SystemGenerator(forcefields=[protein_forcefield,solvation_forcefield],
        barostat=barostat, forcefield_kwargs=forcefield_kwargs, molecules=[ligand],
        small_molecule_forcefield=small_molecule_forcefield)
    
    modeller = app.Modeller(complex_structure.topology, complex_structure.positions)
    modeller.addSolvent(system_generator.forcefield, model='tip3p',
        padding=solvent_padding, ionicStrength=ionic_strength)
    
    system = system_generator.create_system(modeller.topology)
    solvated_structure = parmed.openmm.load_topology(modeller.topology,
        system, xyz=modeller.positions)

    integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    with open(f'{output_prefix}/integrator.xml', 'w') as f:
        f.write(openmm.XmlSerializer.serialize(integrator))


    # minimize and equilibrate
    print('Minimizing...')
    platform = openmm.Platform.getPlatformByName('CUDA')
    platform.setPropertyDefaultValue('CudaDeviceIndex', sys.argv[3])
    context = openmm.Context(system, integrator, platform)
    context.setPositions(modeller.positions)
    openmm.LocalEnergyMinimizer.minimize(context)
    
    print('Equilibrating...')
    integrator.step(nsteps_equil)
    
    state = context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True)
    positions = state.getPositions(asNumpy=True)

    system = context.getSystem()
    system.setDefaultPeriodicBoxVectors(*state.getPeriodicBoxVectors())

    parmed_system = parmed.openmm.load_topology(modeller.topology,
        system, xyz=modeller.positions)
    parmed_system.save(f'{output_prefix}/conf.gro', overwrite=True)
    parmed_system.save(f'{output_prefix}/topol.top', overwrite=True)

  except Exception as e:
    print(e)
    continue
