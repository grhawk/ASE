"""
This module defines an ASE interface to D3H4.
!!!WARNING!!!
This correction can be applied to PM6, DFTB and PM3 and this is defined in the enviromental variable
that compose the command!
"""

import os
import numpy as np
from ase.units import kcal, mol, Angstrom
from ase.calculators.calculator import FileIOCalculator, kpts2mp


class D3H4(FileIOCalculator):
    """ A calculator to compute the D3H4 corrections with ase-FileIOCalculator nomenclature
    """

    if 'D3' in os.environ and 'H4' in os.environ and 'D3FUNC' in os.environ and 'D3DF' in os.environ:
        D3H4_command = os.environ['D3'] + ' struct-d3h4-ase.xyz ' + '-grad -func '+os.environ['D3FUNC']+' -'+os.environ['D3DF']+' > d3h4.out; ' + os.environ['H4'] + ' < struct-d3h4-ase.xyz >> d3h4.out'
    else:
        raise EnvironmentError('1','Some variable is missing: H4, D3, D3FUNC, D3DF')

    command = D3H4_command

    implemented_properties = ['energy','forces']

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='D3H4', atoms=None, kpts=None,
                 **kwargs):
        """Construct a D3H4 calculator.
        """

        self.default_parameters = dict(
            noH4 = False
            )

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        if kpts != None: raise NotImplemented ('D3H4 cannot execute k-points computations')
        if restart != None: raise NotImplemented ('No meaning to restart a D3H4 computation')

        #indexes for the result file
        self.first_time = True
        self.index_energy = None


    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        return system_changes

    def write_input(self, atoms, properties=None, system_changes=None):
        from ase.io import write
        FileIOCalculator.write_input(\
            self, atoms, properties, system_changes)
        write('struct-d3h4-ase.xyz', atoms)

    def read_results(self):
        """ all results are read from label.tag file 
            It will be destroyed after it is read to avoid
            reading it once again after some runtime error """
        from ase.io import read
        from os import remove

        myfile = open('d3h4.out', 'r')
        self.lines = myfile.readlines()
        myfile.close()
        if self.first_time:
            self.first_time = False
            # Energy line index
            for iline, line in enumerate(self.lines):
                d3_estring = 'E6'
                h4_estring = 'Total:'
                if line.find(d3_estring) >= 0:
                    self.index_d3energy = iline
                if line.find(h4_estring) >= 0:
                    self.index_h4energy = iline


        self.read_energy()
        self.read_forces()

            
    def read_energy(self):
        """Read Energy from dDMC tag output file (dDMC.tag)."""
        from ase.units import kcal,mol

        # Energy:
        try:
        # The kcal/mol is to convert everything in eV
            d3_energy = float(self.lines[self.index_d3energy].split()[3]) * kcal/mol
            h4_energy = float(self.lines[self.index_h4energy].split()[1]) * kcal/mol

            # For back compatibility purpose
            if 'H4_correction' in os.environ:
                if os.environ['H4_correction'] == 'no': h4_energy = 0.0
            ################################

            if self.parameters['noH4']: h4_energy = 0.0
            self.results['energy'] = d3_energy + h4_energy
        except:
            raise RuntimeError('Problem in reading energy')

    def read_matrix(self,nfile,string,nats):

        outfile = open(nfile,'r')
        lines = outfile.readlines()
        outfile.close()

        offset = 1

        if string == None: 
            string = lines[0].split()[0]
            offset = 0

        forces = np.zeros((nats, 3), float)

        try:
            for i, line in enumerate(lines):
                if line.find(string) != -1:
                    for j in range(nats):
                        gline = lines[i + j + offset]
                        forces[j,:] = map(float, map(lambda x: x.replace('D','E'), gline.split()))
                    break

            return np.array(forces)

        except:
            raise RuntimeError('Problem in reading forces from '+nfile)

    def read_forces(self):
        """Read Forces from the d3h4.out output file and from dft3_gradient file."""
        from ase.units import kcal, mol, Angstrom
        from ase.units import Hartree, Bohr

        
        nats = len(self.atoms)
        H4_grad = self.read_matrix('d3h4.out','Total gradient', nats) * kcal/mol/Angstrom
        # D3_grad = self.read_matrix('dftd3_gradient',None, nats) * kcal/mol/Angstrom
        D3_grad = self.read_matrix('dftd3_gradient',None, nats) * Hartree / Bohr

        if self.parameters['noH4']: H4_grad = 0.0
        # For back compatibility purpose
        if 'H4_correction' in os.environ:
            if os.environ['H4_correction'] == 'no': H4_grad = 0.0
        ################################
        # print H4_grad
        # print D3_grad
        # print 'ang',kcal/mol

        try:
            self.results['forces'] = D3_grad + H4_grad
        except:
            raise RuntimeError('Problem in reading forces from ')
