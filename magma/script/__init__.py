import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def main():
    "Entry point for magma script"
    magmacommand = MagmaCommand()
    return magmacommand.run()

class MagmaCommand(object):
    """Console script with multiple sub-commands to perform MAGMa calculations"""

    def __init__(self):
        self.parser = argparse.ArgumentParser(description=self.__doc__)
        self.parser.add_argument('--version', action='version', version='%(prog)s '+self.version())
        subparsers = self.parser.add_subparsers(title='Sub-commands')

        sc = subparsers.add_parser("allinone", help=self.allinone.__doc__, description=self.allinone.__doc__)
        # mzxml arguments
        sc.add_argument('mzxml', type=argparse.FileType('r'), help="mzXMl file with MS/MS data")
        sc.add_argument('-p', '--precision', help="precision in Dalton (default: %(default)s)", default=0.01)
        sc.add_argument('-c', '--msfilter', help="cutoff value to filter MS peaks (absolute) (default: %(default)s)", default=2e5)
        sc.add_argument('-d', '--msmsfilter', help="cutoff value to filter MSMS peaks (relative to basepeak) (default: %(default)s)", default=0.1)
        sc.add_argument('-i', '--ionisation', help="Ionisation mode (default: %(default)s)", default="1", choices=["-1", "1"])
        # sygma arguments
        sc.add_argument('seed', type=argparse.FileType('r'), help="File with smiles used as metabolite seeds")
        sc.add_argument('-n', '--nsteps', help="Maximum number of reaction steps (default: %(default)s)", default=2)
        sc.add_argument('-b', '--biotransformations', help="1 and/or 2 for phase 1 and 2 biotransformations (default: %(default)s)", default="12", choices=["1", "2", "12"])
        # output arguments
        sc.add_argument('db', type=argparse.FileType('w'), help="Sqlite database file with results")
        sc.set_defaults(func=self.allinone)

        sc = subparsers.add_parser("metabolize", help=self.metabolize.__doc__, description=self.metabolize.__doc__)
        sc.add_argument('db', type=argparse.FileType('r'), help="Sqlite database file with results")
        # sygma arguments
        sc.add_argument('seed', type=argparse.FileType('r'), help="File with smiles used as metabolite seeds")
        sc.add_argument('-n', '--nsteps', help="Maximum number of reaction steps (default: %(default)s)", default=2)
        sc.add_argument('-b', '--biotransformations', help="1 and/or 2 for phase 1 and 2 biotransformations (default: %(default)s)", default="12", choices=["1", "2", "12"])
        # output arguments
        sc.set_defaults(func=self.metabolize)

        sc = subparsers.add_parser("mzxml", help=self.mzxml.__doc__, description=self.mzxml.__doc__)
        # mzxml arguments
        sc.add_argument('mzxml', type=argparse.FileType('r'), help="mzXMl file with MS/MS data")
        sc.add_argument('-p', '--precision', help="precision in Dalton (default: %(default)s)", default=0.01)
        sc.add_argument('-c', '--msfilter', help="cutoff value to filter MS peaks (absolute) (default: %(default)s)", default=2e5)
        sc.add_argument('-d', '--msmsfilter', help="cutoff value to filter MSMS peaks (relative to basepeak) (default: %(default)s)", default=0.1)
        sc.add_argument('-i', '--ionisation', help="Ionisation mode (default: %(default)s)", default="1", choices=["-1", "1"])
        # output arguments
        sc.add_argument('db', type=argparse.FileType('w'), help="Sqlite database file with results")
        sc.set_defaults(func=self.mzxml)

        sc = subparsers.add_parser("sd2smiles", help=self.sd2smiles.__doc__, description=self.sd2smiles.__doc__)
        sc.add_argument('input', type=argparse.FileType('r'), help="Sd file")
        sc.add_argument('output', type=argparse.FileType('w'), help="File with smiles which can be used as metabolite seeds")
        sc.set_defaults(func=self.sd2smiles)

        sc = subparsers.add_parser("smiles2sd", help=self.smiles2sd.__doc__, description=self.smiles2sd.__doc__)
        sc.add_argument('input', type=argparse.FileType('r'), help="File with smiles which can be used as metabolite seeds")
        sc.add_argument('output', type=argparse.FileType('w'), help="Sd file")
        sc.add_argument('-d', '--depiction', help="Depiction (default: %(default)s)", default="1D", choices=["1D", "2D", "3D"])
        sc.set_defaults(func=self.smiles2sd)

    def version(self):
        return '1.0' # TODO move to main magma package and reuse in setup.py so version is specified in one place

    def allinone(self, args):
        """Reads metabolite seed file and MS/MS datafile, generates metabolites and matches them to peaks"""
        print 'allinone'
        print args
        # TODO implement

    def metabolize(self, args):
        """Reads metabolite seed file and existing result database, generates metabolites and matches them to peaks"""
        print 'metabolize'
        print args
        # TODO implement

    def mzxml(self, args):
        """Reads MS/MS datafile"""
        print 'mzxml'
        print args
        # TODO implement

    def sd2smiles(self, args):
        """ Convert sd file to smiles """
        mols = Chem.SDMolSupplier(args.input)
        for mol in mols:
            args.output.writeln(
                               '%s|%s'.format(
                                              Chem.MolToSmiles(mol),
                                              mol.GetProp('_Name')
                                              )
                               )

    def smiles2sd(self, args):
        """ Convert smiles to sd file"""
        w = Chem.SDWriter(args.output)
        for line in args.input:
            line = line.strip()
            (smilestring, molname) = line.split('|')
            mol = Chem.MolFromSmiles(smilestring)
            if (args.depiction == '2D'):
                AllChem.Compute2DCoords(mol)
            elif (args.depiction == '3D'):
                AllChem.EmbedMolecule(mol)
                AllChem.UFFOptimizeMolecule(mol)
            mol.SetProp('_Name', molname)
            w.write(mol)

    def run(self):
        """Parse arguments"""
        args =  self.parser.parse_args()
        return args.func(args)
