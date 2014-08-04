import os

import numpy as np

#from ase.calculators.calculator import FileIOCalculator, kpts2mp
from ase.calculators.dftb import Dftb
from ase.calculators.ddmc import dDMC

from ase.calculators.calculator import FileIOCalculator, kpts2mp

class DftbdDMC(FileIOCalculator):
    """ A dftb+-dDMC calculator with ase-FileIOCalculator nomenclature
    """

    implemented_properties = ['energy', 'forces']#, 'charges']

    def __init__(self, label='dftbddmc', atoms=None, dftbdict={}, ddmcdict={}, **kwargs):
        """Construct a DFTB+ and a dDMC calculator.
        """

        os.environ['ASE_DFTBDDMC_COMMAND'] = ''

        self.label = label
        self.dftbdict = dftbdict
        self.ddmcdict = ddmcdict
        
        ddmc_default_parameters = {}
        ddmc_default_parameters['param_a'] = 1.85705835084132
        ddmc_default_parameters['param_b'] = 1.01824853175310
        ddmc_default_parameters['param_c'] = 23.0
        ddmc_default_parameters['dftype'] = 3
        ddmc_default_parameters['tagtype'] = 'dftbp'

        for param in ddmc_default_parameters.keys():
            if not param in self.ddmcdict.keys():
                self.ddmcdict[param] = ddmc_default_parameters[param]


        restart = None
        ignore_bad_restart_file = None
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        self.dftb_calc = Dftb(label=self.label)
        self.ddmc_calc = dDMC(label=self.label)


    def write_input(self, atoms, properties=None, systems_changes=None):
        self.dftb_calc.parameters.update(self.dftbdict)
        self.ddmc_calc.parameters.update(self.ddmcdict)

        
    def read_results(self):
        self.read_energy()
        self.read_forces()
        
    def read_energy(self):
        results=[]
        """Read Energy from dftb output file (results.tag)."""
        from ase.units import Hartree
        # Energy:
        try:
            self.atoms.set_calculator(self.dftb_calc)
            dftb_energy = self.atoms.get_potential_energy()
            self.atoms.set_calculator(self.ddmc_calc)
            ddmc_energy = self.atoms.get_potential_energy()
            self.results['energy'] = ddmc_energy+dftb_energy
        except:
            raise RuntimeError('Problem in reading energy')

    def read_forces(self):
        """Read Forces from dftb output file (results.tag)."""
        from ase.units import Hartree, Bohr

        try:
            self.atoms.set_calculator(self.dftb_calc)
            dftb_forces = self.atoms.get_forces()
            self.atoms.set_calculator(self.ddmc_calc)
            ddmc_forces = self.atoms.get_forces()
            
            self.results['forces'] = dftb_forces+ddmc_forces

        except:
            raise RuntimeError('Problem in reading forces')
