import os

import numpy as np

#from ase.calculators.calculator import FileIOCalculator, kpts2mp
from ase.calculators.mopac import Mopac
from ase.calculators.ddmc import dDMC

from ase.calculators.calculator import FileIOCalculator, kpts2mp

class MopacdDMC(FileIOCalculator):
    """ A MOPAC-dDMC calculator with ase-FileIOCalculator nomenclature
    """

    implemented_properties = ['energy', 'forces']#, 'charges']

    def __init__(self, label='mopacddmc', atoms=None, mopacdict={}, ddmcdict={}, **kwargs):
        """Construct a MOPAC and a dDMC calculator.
        """

        os.environ['ASE_MOPACDDMC_COMMAND'] = ''

        self.label = label
        self.mopacdict = mopacdict
        self.ddmcdict = ddmcdict
        
        self.ddmcdict['param_a'] = 1.5880952
        self.ddmcdict['param_b'] = 0.2533719
        self.ddmcdict['param_c'] = 23.0
        self.ddmcdict['dftype'] = 3
        self.ddmcdict['tagtype'] = 'column'


        restart = None
        ignore_bad_restart_file = None
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)


    def write_input(self, atoms, properties=None, systems_changes=None):
        self.gaus_calc = Gaussian(label=self.label)
        self.gaus_calc.parameters.update(self.gausdict)
        self.ddmc_calc = dDMC(label=self.label)
        self.ddmc_calc.parameters.update(self.ddmcdict)

        
    def read_results(self):
        self.read_energy()
        self.read_forces()
#        self.read_charges()
        
    def read_energy(self):
        results=[]
        """Read Energy from dftb output file (results.tag)."""
        from ase.units import Hartree
        # Energy:
        try:
            self.atoms.set_calculator(self.mopac_calc)
            mopac_energy = self.atoms.get_potential_energy()
            mopac_pop = self.atoms.get_charges()
            tag = open(self.label+'.tag','w')
            tag.write(str(len(self.atoms))+'\n')
            for p in mopac_pop:
                tag.write(str(p)+'\n')
            tag.close()
            self.atoms.set_calculator(self.ddmc_calc)
            ddmc_energy = self.atoms.get_potential_energy()
            self.results['energy'] = ddmc_energy+mopac_energy
        except:
            raise RuntimeError('Problem in reading energy')

    def read_forces(self):
        """Read Forces from dftb output file (results.tag)."""
        from ase.units import Hartree, Bohr

        try:
            self.atoms.set_calculator(self.mopac_calc)
            mopac_forces = self.atoms.get_forces()
            self.atoms.set_calculator(self.ddmc_calc)
            ddmc_forces = self.atoms.get_forces()
            
            self.results['forces'] = mopac_forces+ddmc_forces

        except:
            raise RuntimeError('Problem in reading forces')
        
    def read_charges(self):
        try:
            self.atoms.set_calculator(self.dftb_calc)
            self.results['charges'] = self.atoms.get_charges()



        except:
            raise RuntimeError('Problem in reading charges')

