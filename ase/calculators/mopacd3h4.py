import os

import numpy as np

#from ase.calculators.calculator import FileIOCalculator, kpts2mp
from ase.calculators.mopac import Mopac
from ase.calculators.d3h4 import D3H4

from ase.calculators.calculator import FileIOCalculator, kpts2mp

class MopacD3H4(FileIOCalculator):
    """ A MOPAC-dDMC calculator with ase-FileIOCalculator nomenclature
    """

    implemented_properties = ['energy', 'forces']#, 'charges']

    def __init__(self, label='mopacddmc', atoms=None, mopacdict={}, d3h4dict={}, **kwargs):
        """Construct a MOPAC and a dDMC calculator.
        """

        os.environ['ASE_MOPACD3H4_COMMAND'] = ''

        self.label = label
        self.mopacdict = mopacdict
        self.d3h4dict = d3h4dict
        
        restart = None
        ignore_bad_restart_file = None
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)


    def write_input(self, atoms, properties=None, systems_changes=None):
        self.mopac_calc = Mopac(label=self.label)
        self.mopac_calc.parameters.update(self.mopacdict)
        self.d3h4_calc = D3H4(label=self.label)
        self.d3h4_calc.parameters.update(self.d3h4dict)

        
    def read_results(self):
        self.read_energy()
        self.read_forces()
        
    def read_energy(self):
        results=[]
        """Read Energy from dftb output file (results.tag)."""
        from ase.units import Hartree
        # Energy:
        try:
            self.atoms.set_calculator(self.mopac_calc)
            mopac_energy = self.atoms.get_potential_energy()
            self.atoms.set_calculator(self.d3h4_calc)
            d3h4_energy = self.atoms.get_potential_energy()
            self.results['energy'] = d3h4_energy+mopac_energy
        except:
            raise RuntimeError('Problem in reading energy')

    def read_forces(self):
        """Read Forces from dftb output file (results.tag)."""
        from ase.units import Hartree, Bohr

        try:
            self.atoms.set_calculator(self.mopac_calc)
            mopac_forces = self.atoms.get_forces()
            self.atoms.set_calculator(self.d3h4_calc)
            d3h4_forces = self.atoms.get_forces()
            
            self.results['forces'] = mopac_forces+d3h4_forces

        except:
            raise RuntimeError('Problem in reading forces')
        
