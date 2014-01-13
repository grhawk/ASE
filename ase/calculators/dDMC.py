"""This module defines an ASE interface to DftbPlus

http://http://www.dftb-plus.info//
http://www.dftb.org/

markus.kaukonen@iki.fi

The file 'geom.out.gen' contains the input and output geometry
and it will be updated during the dftb calculations.

If restart == None
                   it is assumed that a new input file 'dftb_hsd.in'
                   will be written by ase using default keywords
                   and the ones given by the user.

If restart != None
                   it is assumed that keywords are in file restart

The keywords are given, for instance, as follows::

    Hamiltonian_SCC ='YES',
    Hamiltonian_SCCTolerance = 1.0E-008,
    Hamiltonian_MaxAngularMomentum = '',
    Hamiltonian_MaxAngularMomentum_O = '"p"',
    Hamiltonian_MaxAngularMomentum_H = '"s"',
    Hamiltonian_InitialCharges_ = '',
    Hamiltonian_InitialCharges_AllAtomCharges_ = '',
    Hamiltonian_InitialCharges_AllAtomCharges_1 = -0.88081627,
    Hamiltonian_InitialCharges_AllAtomCharges_2 = 0.44040813,
    Hamiltonian_InitialCharges_AllAtomCharges_3 = 0.44040813,

"""

import os

import numpy as np

from ase.calculators.calculator import FileIOCalculator, kpts2mp


class dDMC(FileIOCalculator):
    """ A calculator to compute the dDMC corrections with ase-FileIOCalculator nomenclature
    """
    if 'DFTB_COMMAND' in os.environ:
        dftb_command = os.environ['DFTB_COMMAND'] + ' > DFTB.out'
    else:
        dftb_command = 'dftb+ > PREFIX.out'

    if 'dDMC_COMMAND' in os.environ:
        dDMC_command = os.environ['DFTB_COMMAND'] + ' > dDMC.out'
    else:
        raise MissingEnviromentVariable('DFTB_COMMAND')

    command = dDMC_command

    implemented_properties = ['energy']

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='dDMC', atoms=None, kpts=None,
                 **kwargs):
        """Construct a dDMC calculator.
        """

        self.default_parameters = dict(
            tag = label+'.tag',
            debugflag = 'down'
            )

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        if kpts != None: raise NotImplemented ('dDMC cannot execute k-points computations')
        if restart != None: raise NotImplemented ('No meaning to restart a dDMC computation')

        #indexes for the result file
        self.first_time = True
        self.index_energy = None

    def write_dDMC_in(self):
        """ Write the innput file for the dDMC calculation.
            Geometry is taken always from the file 'label.xyz'.
        """

        outfile = open('ddmc.in', 'w')
        outfile.write('# This file is prepared by python-ASE')

        #--------MAIN KEYWORDS-------
        for key,value in sorted(self.parameters.items()):
            outfile.write(key+' = '+parameters[key])
        

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()
            self.write_dDMC_in()

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        return system_changes

    def write_input(self, atoms, properties=None, system_changes=None):
        from ase.io import write
        FileIOCalculator.write_input(\
            self, atoms, properties, system_changes)
        self.write_dDMC_in()
        write(label+'.xyz', atoms)

    def read_results(self):
        """ all results are read from results.tag file 
            It will be destroyed after it is read to avoid
            reading it once again after some runtime error """
        from ase.io import read
        from os import remove

        myfile = open('results.tag', 'r')
        self.lines = myfile.readlines()
        myfile.close()
        if self.first_time:
            self.first_time = False
            # Energy line index
            for iline, line in enumerate(self.lines):
                estring = 'total_energy'
                if line.find(estring) >= 0:
                    self.index_energy = iline + 1
                    break
            # Force line indexes
            for iline, line in enumerate(self.lines):
                fstring = 'forces   '
                if line.find(fstring) >= 0:
                    self.index_force_begin = iline + 1
                    line1 = line.replace(':', ',')
                    self.index_force_end = iline + 1 + \
                        int(line1.split(',')[-1])
                    break
            # Charge line indexes
            for iline, line in enumerate(self.lines):
                fstring = 'net_atomic_charges'
                if line.find(fstring) >= 0:
                    self.index_charge_begin = iline + 1
                    line1 = line.replace(':', ',')
                    natoms = int(line1.split(',')[-1])
                    mod = natoms%3
                    self.index_charge_end = iline + 1 + \
                         natoms/3 + mod
                    break


        self.read_energy()
        self.read_charges()
        # read geometry from file in case dftb+ has done steps
        # to move atoms, in that case forces are not read
        if int(self.parameters['Driver_MaxSteps']) > 0:
            self.atoms = read('geo_end.gen')
            self.results['forces'] = np.zeros([len(self.state), 3])
        else:
            self.read_forces()
        os.remove('results.tag')
            
    def read_energy(self):
        """Read Energy from dftb output file (results.tag)."""
        from ase.units import Hartree

        # Energy:
        try:
            energy = float(self.lines[self.index_energy].split()[0]) * Hartree
            self.results['energy'] = energy
        except:
            raise RuntimeError('Problem in reading energy')

    def read_forces(self):
        """Read Forces from dftb output file (results.tag)."""
        from ase.units import Hartree, Bohr

        try:
            gradients = []
            for j in range(self.index_force_begin, self.index_force_end):
                word = self.lines[j].split()
                gradients.append([float(word[k]) for k in range(0, 3)])

            self.results['forces'] = np.array(gradients) * Hartree / Bohr

        except:
            raise RuntimeError('Problem in reading forces')
        
    def read_charges(self):
        try:
            charges = []
            for j in range(self.index_charge_begin, self.index_charge_end):
                word = self.lines[j].split()
                charges.append([float(word[k]) for k in range(len(word))])

            self.results['charges'] = np.array(charges)

        except:
            raise RuntimeError('Problem in reading charges')

