"""This module defines an ASE interface to dDMC
"""

import os

import numpy as np

from ase.calculators.calculator import FileIOCalculator, kpts2mp


class dDMC(FileIOCalculator):
    """ A calculator to compute the dDMC corrections with ase-FileIOCalculator nomenclature
    """
    if 'dDMC_COMMAND' in os.environ:
        dDMC_command = os.environ['dDMC_COMMAND'] + ' < '+ 'ddmc.in' +' > dDMC.out'
    else:
        raise EnvironmentError('1','dDMC_COMMAND has to be defined')

    command = dDMC_command

    implemented_properties = ['energy','forces']

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='dDMC', atoms=None, kpts=None,
                 **kwargs):
        """Construct a dDMC calculator.
        """

        self.default_parameters = dict(
            tag = label+'.tag',
            debugflag = 'down',
            geometry = label+'.xyz',
            atomdata = 'notImportant',
            tagtype = 'dftbp',
            gradient = 'UP'
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
        outfile.write('# This file is prepared by python-ASE\n')
        parfile_open = False

        #--------MAIN KEYWORDS-------
        for key,value in sorted(self.parameters.items()):
            if key.find('param_') >=0:

                if parfile_open == False:
                    parfile = open('parameters.dat', 'w')
                    parfile.write('# This file is prepared by python-ASE\n')
                    parfile_open = True
                    
                parfile.write('%12.6f\n' % float(value))
            else:
                outfile.write(str(key)+' = '+str(value)+'\n')
        

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
        write(self.label+'.xyz', atoms)

    def read_results(self):
        """ all results are read from label.tag file 
            It will be destroyed after it is read to avoid
            reading it once again after some runtime error """
        from ase.io import read
        from os import remove

        myfile = open('dDMC.tag', 'r')
        self.lines = myfile.readlines()
        myfile.close()
        if self.first_time:
            self.first_time = False
            # Energy line index
            for iline, line in enumerate(self.lines):
                estring = 'correction_energy'
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
            # # Charge line indexes
            # for iline, line in enumerate(self.lines):
            #     fstring = 'net_atomic_charges'
            #     if line.find(fstring) >= 0:
            #         self.index_charge_begin = iline + 1
            #         line1 = line.replace(':', ',')
            #         natoms = int(line1.split(',')[-1])
            #         mod = natoms%3
            #         self.index_charge_end = iline + 1 + \
            #              natoms/3 + mod
            #         break


        self.read_energy()
        self.read_forces()
        # # read geometry from file in case dftb+ has done steps
        # # to move atoms, in that case forces are not read
        # if int(self.parameters['Driver_MaxSteps']) > 0:
        #     self.atoms = read('geo_end.gen')
        #     self.results['forces'] = np.zeros([len(self.state), 3])
        # else:
        #     self.read_forces()
#        os.remove('dDMC.tag')
            
    def read_energy(self):
        """Read Energy from dDMC tag output file (dDMC.tag)."""
        from ase.units import Hartree

        # Energy:
        try:
            energy = float(self.lines[self.index_energy]) * Hartree
            self.results['energy'] = energy
        except:
            raise RuntimeError('Problem in reading energy')

    def read_forces(self):
        """Read Forces from dDMC output file (results.tag)."""
        from ase.units import Hartree, Bohr

        try:
            gradients = []
            for j in range(self.index_force_begin, self.index_force_end):
                word = self.lines[j].split()
                gradients.append([float(word[k]) for k in range(0, 3)])
                
            self.results['forces'] = np.array(gradients) * Hartree / Bohr

        except:
            raise RuntimeError('Problem in reading forces')
        
    # def read_charges(self):
    #     try:
    #         charges = []
    #         for j in range(self.index_charge_begin, self.index_charge_end):
    #             word = self.lines[j].split()
    #             charges.append([float(word[k]) for k in range(len(word))])

    #         self.results['charges'] = np.array(charges)

    #     except:
    #         raise RuntimeError('Problem in reading charges')

