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

        ddmc_default_parameters = {}
        
        ddmc_default_parameters['param_a'] = 1.53014262236515
        ddmc_default_parameters['param_b'] = 1.04216259309481
        ddmc_default_parameters['param_c'] = 23.0
        ddmc_default_parameters['param_d'] = 1.84166256639050
        ddmc_default_parameters['dftype'] = 4
        ddmc_default_parameters['tagtype'] = 'column'

        for param in ddmc_default_paramters.keys():
            if not param in self.ddmcdict.keys():
                self.ddmcdict[param] = ddmc_default_parameters[param]


        restart = None
        ignore_bad_restart_file = None
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)


    def write_input(self, atoms, properties=None, systems_changes=None):
        from ase.io import write
        atoms.write(self.label+'.xyz','xyz')
        self.mopac_calc = Mopac(label=self.label)
        self.mopac_calc.parameters.update(self.mopacdict)
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

