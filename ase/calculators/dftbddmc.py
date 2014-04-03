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

        # from ase.dft.kpoints import monkhorst_pack

        # if 'DFTB_PREFIX' in os.environ:
        #     slako_dir = os.environ['DFTB_PREFIX']
        # else:
        #     slako_dir = './'

        # self.default_parameters = dict(
        #     Hamiltonian_='DFTB',
        #     Driver_='ConjugateGradient',
        #     Driver_MaxForceComponent='1E-4',
        #     Driver_MaxSteps=0,
        #     Hamiltonian_SlaterKosterFiles_='Type2FileNames',
        #     Hamiltonian_SlaterKosterFiles_Prefix=slako_dir,
        #     Hamiltonian_SlaterKosterFiles_Separator='"-"',
        #     Hamiltonian_SlaterKosterFiles_Suffix='".skf"',
        #     Hamiltonian_SCC = 'No',
        #     )

        restart = None
        ignore_bad_restart_file = None
        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        # self.kpts = kpts
        # # kpoint stuff by ase
        # if self.kpts != None:
        #     mpgrid = kpts2mp(atoms, self.kpts)
        #     mp = monkhorst_pack(mpgrid)
        #     initkey = 'Hamiltonian_KPointsAndWeights'
        #     self.parameters[initkey + '_'] = ''
        #     for i, imp in enumerate(mp):
        #         key = initkey + '_empty' + str(i)
        #         self.parameters[key] = str(mp[i]).strip('[]') + ' 1.0'

        #the input file written only once
        # if restart == None:
        #     self.write_dftb_in()
        # else:
        #     if os.path.exists(restart):
        #         os.system('cp ' + restart + ' dftb_in.hsd')
        #     if not os.path.exists('dftb_in.hsd'):
        #         raise IOError('No file "dftb_in.hsd", use restart=None')

        #indexes for the result file
        # self.first_time = True
        # self.index_energy = None
        # self.index_force_begin = None
        # self.index_force_end = None
        # self.index_charge_begin = None
        # self.index_charge_end = None

        # self.dftb_calc = Dftb(label=label)
        # self.dftb_calc.parameters.update(dftbdict)
        # self.dftb_energy = atoms.get_potential_energy()
        # self.ddmc_calc = dDMC(atoms=atoms,label=label)
        # self.ddmc_calc.parameters.update(ddmcdict)
        # self.atoms = atoms


    def write_input(self, atoms, properties=None, systems_changes=None):
        self.dftb_calc = Dftb(label=self.label)
        self.dftb_calc.parameters.update(self.dftbdict)
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
        
    def read_charges(self):
        try:
            self.atoms.set_calculator(self.dftb_calc)
            self.results['charges'] = self.atoms.get_charges()



        except:
            raise RuntimeError('Problem in reading charges')

