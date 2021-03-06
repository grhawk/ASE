"""
Gaussian calculator for ASE written by:

    Glen R. Jenness
    University of Wisconsin - Madison

Based off of code written by:

    Glen R. Jenness
    Kuang Yu
    Torsten Kerber, Ecole normale superieure de Lyon (*)
    Paul Fleurat-Lessard, Ecole normale superieure de Lyon (*)
    Martin Krupicka

(*) This work is supported by Award No. UK-C0017, made by King Abdullah
University of Science and Technology (KAUST), Saudi Arabia.

See accompanying license files for details.
"""
import os

from ase.calculators.calculator import FileIOCalculator

"""
Gaussian has two generic classes of keywords:  link0 and route.
Since both types of keywords have different input styles, we will
distinguish between both types, dividing each type into str's, int's
etc.

For more information on the Link0 commands see:
    http://www.gaussian.com/g_tech/g_ur/k_link0.htm
For more information on the route section keywords, see:
    http://www.gaussian.com/g_tech/g_ur/l_keywords09.htm
"""
link0_keys = [\
              'chk',
              'mem',
              'rwf',
              'int',
              'd2e',
              'lindaworkers',
              'kjob',
              'subst',
              'save',
              'nosave',
              'nprocshared',
              'nproc',
             ]

# This one is a little strange.  Gaussian has several keywords where you just
# specify the keyword, but the keyword itself has several options.
# Ex:  Opt, Opt=QST2, Opt=Conical, etc.
# These keywords are given here.
route_self_keys = ['opt',
                   'force',
                   'freq',
                   'complex',
                   'fmm',
                   'genchk',
                   'polar',
                   'prop',
                   'pseudo',
                   'restart',
                   'scan',
                   'scrf',
                   'sp',
                   'sparse',
                   'stable',
                   'volume',
                  ]

route_keys = [\
# int keys
# Multiplicity and charge are not really route keywords, but we will
# put them here anyways
              'cachesize',
              'cbsextrapolate',
              'constants',
# str keys
              'functional',
              'maxdisk',
              'cphf',
              'density',
              'densityfit',
              'ept',
              'field',
              'geom',
              'guess',
              'gvb',
              'integral',
              'irc',
              'ircmax',
              'name',
              'nmr',
              'nodensityfit',
              'oniom',
              'output',
              'punch',
              'scf',
              'symmetry',
              'td',
              'units',
# Float keys
              'pressure',
              'scale',
              'temperature',
             ]


class Gaussian(FileIOCalculator):
    """
    Gaussian calculator
    """
    name = 'Gaussian'

    implemented_properties = ['energy', 'forces', 'charges', 'dipole']
    command = 'g09 < PREFIX.com > PREFIX.log'

    def __init__(self, restart=None, ignore_bad_restart_file=False,
                 label='g09', atoms=None, scratch=None, ioplist=list(),
                 basisfile=None, **kwargs):

        """Constructs a Gaussian-calculator object.

        """

        FileIOCalculator.__init__(self, restart, ignore_bad_restart_file,
                                  label, atoms, **kwargs)

        self.ioplist = ioplist
        self.scratch = scratch
        self.basisfile = basisfile

    def set(self, **kwargs):
        changed_parameters = FileIOCalculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()
        return changed_parameters

    def initialize(self, atoms):
# Set some default behavior
        if ('multiplicity' not in self.parameters):
            self.parameters['multiplicity'] = 1

        if ('charge' not in self.parameters):
            self.parameters['charge'] = 0

        if ('method' not in self.parameters):
            self.parameters['method'] = 'hf'

        if ('basis' not in self.parameters):
            self.parameters['basis'] = '6-31g*'

        if ('force' not in self.parameters):
            self.parameters['force'] = 'force'

        self.converged = None

    def write_input(self, atoms, properties=None, system_changes=None):
        """Writes the input file"""
        FileIOCalculator.write_input(self, atoms, properties, system_changes)
        self.initialize(atoms)

        filename = self.label + '.com'
        inputfile = open(filename, 'w')

        link0 = str()
        route = '#p %s/%s' % (self.parameters['method'],
                              self.parameters['basis'])

        for key, val in self.parameters.items():
            if key.lower() in link0_keys:
                link0 += ('%%%s=%s\n' % (key, val))
            elif key.lower() in route_self_keys:
                if (val.lower() == key.lower()):
                    route += (' ' + val)
                else:
                    if ',' in val:
                        route += ' %s(%s)' % (key, val)
                    else:
                        route += ' %s=%s' % (key, val)
            elif key.lower() in route_keys:
                route += ' %s=%s' % (key, val)

        if (self.ioplist):
            route += ' IOp('
            for iop in self.ioplist:
                route += (' ' + iop)
                if (len(self.ioplist) > 1) and (iop != len(self.ioplist) - 1):
                    route += ','
            route += ')'

        inputfile.write(link0)
        inputfile.write(route)
        inputfile.write(' \n\n')
        inputfile.write('Gaussian input prepared by ASE\n\n')
        inputfile.write('%i %i\n' % (self.parameters['charge'],
                                     self.parameters['multiplicity']))

        symbols = atoms.get_chemical_symbols()
        coordinates = atoms.get_positions()
        for i in range(len(atoms)):
            inputfile.write('%-10s' % symbols[i])
            for j in range(3):
                inputfile.write('%20.10f' % coordinates[i, j])
            inputfile.write('\n')

        inputfile.write('\n')

        if ('gen' in self.parameters['basis'].lower()):
            if (self.basisfile is None):
                raise RuntimeError('Please set basisfile.')
            elif (not os.path.isfile(self.basisfile)):
                raise RuntimeError('Basis file %s does not exist.' \
                % self.basisfile)
            else:
                f2 = open(self.basisfile, 'r')
                inputfile.write(f2.read())
                f2.close()

        if atoms.get_pbc().any():
            cell = atoms.get_cell()
            line = str()
            for v in cell:
                line += 'TV %20.10f%20.10f%20.10f\n' % (v[0], v[1], v[2])
            inputfile.write(line)

        inputfile.write('\n\n')

        inputfile.close()

    def read(self, label):
        """Used to read the results of a previous calculation if restarting"""
        FileIOCalculator.read(self, label)

    def read_results(self):
        """Reads the output file using GaussianReader"""
        from ase.io.gaussian import read_gaussian_out
        filename = self.label + '.log'

        self.results['energy'] = read_gaussian_out(filename, quantity='energy')
        self.results['forces'] = read_gaussian_out(filename, quantity='forces')
        self.results['dipole'] = read_gaussian_out(filename, quantity='dipole')

    def clean(self):
        """Cleans up from a previous run"""
        extensions = ['.chk', '.com', '.log']

        for ext in extensions:
            f = self.label + ext
            try:
                if (self.directory is not None):
                    os.remove(os.path.join(self.directory, f))
                else:
                    os.remove(f)
            except OSError:
                pass

    def read_convergence(self):
        """Determines if calculations converged"""
        converged = False

        gauss_dir = os.environ['GAUSS_EXEDIR']
        test = '(Enter ' + gauss_dir + '/l9999.exe)'

        f = open(self.label + '.log', 'r')
        lines = f.readlines()
        f.close()

        for line in lines:
            if (line.rfind(test) > -1):
                converged = True
            else:
                converged = False

        return converged

    def get_version(self):
        return self.read_output(self.label + '.log', 'version')
