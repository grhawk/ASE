"""
Version 2012/08/20, Torsten Kerber

Contributors:
  Torsten Kerber, Ecole normale superieure de Lyon:
  Paul Fleurat-Lessard, Ecole normale superieure de Lyon
  based on a script by Rosa Bulo, Ecole normale superieure de Lyon

This work is supported by Award No. UK-C0017, made by King Abdullah
University of Science and Technology (KAUST), Saudi Arabia

See accompanying license files for details.
"""
import os
import numpy as np

from ase.units import kcal, mol
from ase.calculators.calculator import FileIOCalculator

class Mopac(FileIOCalculator):
    name = 'MOPAC'

    """Command to run mopac. """

    command = os.environ['MOPAC_COMMAND']+' mopac.mop 2> .mopaclog'

    implemented_properties = ['energy', 'forces', 'charges']

    def __init__(self, label='mopac', restart=None, atoms=None, ignore_bad_restart_file=None,**kwargs):
        
        # set initial values
        self.default_parameters = dict(
            restart = 0,
            spin = 1,
            opt = False,
            functional = 'PM6',
            job_type = ['NOANCI','GRADIENTS', '1SCF'],
            relscf = 0.1,
            charge = 0
        )
        
        # save label and atoms
        self.label = label
        self.atoms = atoms

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file, label, atoms, **kwargs)

        #the input file written only once
        # self.write_mopac_in(self.atoms)

        # initialize the results
        self.version = None
        self.energy_zero = None
        self.energy_free = None
        self.forces = None
        self.charges = None
        self.stress = None
        
        # initialize the results
        self.occupations = None
        
        
    def set(self, **kwargs):
        """
        Sets the parameters on the according keywords
        Raises RuntimeError when wrong keyword is provided
        """
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()
#            self.write_input(self.atoms)

    def get_version(self):
        return self.version

    def write_input(self,atoms,properties=None,system_changes=None):
        """
        Writes the files that have to be written each timestep
        """

        FileIOCalculator.write_input(self,atoms,properties,system_changes)
        
        # start the input
        mopac_input = ''

        #write functional and job_type
        mopac_input += self.parameters['functional']+' '
        for value in self.parameters['job_type']:
            mopac_input += value + ' '
        
        mopac_input += 'RELSCF=' + str(self.parameters['relscf']) + ' '
            
        #write charge
        charge = sum(atoms.get_initial_charges())
        charge = max([charge,self.parameters['charge']])

        if charge != 0:
            mopac_input += 'CHARGE=%i ' % (int(charge))
        
        #write spin
        spin = int(self.parameters['spin'])
        spin_keywords = {
            1: 'SINGLET',
            2: 'DOUBLET',
            3: 'TRIPLET',
        }
        mopac_input += spin_keywords[spin]+' '

        #input down
        mopac_input += '\n'
        mopac_input += 'Title: ASE job\n\n'

        symbols = atoms.get_chemical_symbols()
        coordinates = atoms.get_positions()
        for i in range(len(atoms)):
            mopac_input += ('%-10s' % symbols[i])
            for j in range(3):
                mopac_input += ('%20.10f' % coordinates[i, j])
            mopac_input += ('\n')

        if atoms.pbc.any():
            for v in atoms.get_cell():
                mopac_input += 'Tv %8.3f %8.3f %8.3f\n' % (v[0], v[1], v[2])
        
        # write input
        myfile = open('mopac.mop', 'w')
        myfile.write(mopac_input)
        myfile.close()

    def read_results(self):
        """
        Writes input in label.mop
        Runs MOPAC
        Reads Version, Energy and Forces
        """
        # set the input and output file name
        finput = self.label + '.mop'
        foutput = self.label+'.out'

        self.version = self.read_version()

        energy = self.read_energy()
        self.energy_zero = energy
        self.energy_free = energy
        
        self.forces = self.read_forces()
        self.charges = self.read_charges()

        os.rename('mopac.mop',finput)
        os.rename('mopac.out',foutput)

    def read_version(self):
        """
        Reads the MOPAC version string from the second line
        """
        version = 'unknown'
        lines = open('mopac.out').readlines()
        for line in lines:
            if "  Version" in line:
                version = line.split()[-2]
                break
        return version

    def read_energy(self):
        """
        Reads the ENERGY from the output file (HEAT of FORMATION in kcal / mol)
        Raises RuntimeError if no energy was found
        """
        outfile = open('mopac.out')
        lines = outfile.readlines()
        outfile.close()

        energy = None
        for line in lines:
            # if line.find('TOTAL ENERGY') != -1:
            #     words = line.split()
            #     energy = words[3]
            if line.find('HEAT OF FORMATION') != -1:
                words = line.split()
                energy = words[5]
            # if line.find('H.o.F. per unit cell') != -1:
            #     words = line.split()
            #     energy = words[5]
            if line.find('UNABLE TO ACHIEVE SELF-CONSISTENCE') != -1:
                energy = None
        if energy is None:
            raise RuntimeError('MOPAC: could not find total energy')
        
        try:
            energy = float(energy)
            energy *= (kcal / mol)
            self.results['energy'] = energy
        except:
            raise RuntimeError('Problem in reading energy')
    

    def read_forces(self):
        """
        Reads the FORCES from the output file
        """
        outfile = open('mopac.out')
        lines = outfile.readlines()
        outfile.close()

        nats = len(self.atoms)
        forces = np.zeros((nats, 3), float)

        try:
            for i, line in enumerate(lines):
                if line.find('GRADIENT\n') != -1:
                    for j in range(nats * 3):
                        gline = lines[i + j + 1]
                        forces[j / 3, j % 3] = float(gline[49:62])
                    break

            self.results['forces'] = np.array(forces) * kcal/mol

        except:
            raise RuntimeError('Problem in reading forces')


    def read_charges(self):
        """
        Reads the CHARGES from the output file
        """
        outfile = open('mopac.out')
        lines = outfile.readlines()
        outfile.close()

        nats = len(self.atoms)
        charges = np.zeros(nats, float)
        
        try:
            for i, line in enumerate(lines):
                if line.find('ATOM NO.   TYPE          CHARGE      No. of ELECS.   s-Pop       p-Pop') != -1:
                    for j in range(nats):
                        gline = lines[i + j + 1]
                        charges[j] = float(gline[37:50])
                    break

            self.results['charges'] = np.array(charges)
            
        except:
            raise RuntimeError('Problem in reading charges')

        
    def check_state(self, atoms):
        system_changes = FileIOCalculator.check_state(self, atoms)
        return system_changes

