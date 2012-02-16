import argparse
from rdkit import Chem
from rdkit.Chem import AllChem
import magma

def main():
    "Entry point for magma script"
    magmacommand = MagmaCommand()
    return magmacommand.run()

class MagmaCommand(object):
    """Console script with multiple sub-commands to perform MAGMa calculations"""

    def __init__(self):
        self.parser = argparse.ArgumentParser(description=self.__doc__)
        self.parser.add_argument('--version', action='version', version='%(prog)s ' + self.version())
        subparsers = self.parser.add_subparsers(title='Sub-commands')

        sc = subparsers.add_parser("allinone", help=self.allinone.__doc__, description=self.allinone.__doc__)
        # mzxml arguments
        sc.add_argument('mzxml', type=argparse.FileType('r'), help="mzXMl file with MS/MS data")
        sc.add_argument('-p', '--mz_precision', help="precision in Dalton (default: %(default)s)", default=0.01,type=float)
        sc.add_argument('-c', '--ms_intensity_cutoff', help="cutoff value to filter MS peaks (absolute) (default: %(default)s)", default=2e5,type=float)
        sc.add_argument('-d', '--msms_intensity_cutoff', help="cutoff value to filter MSMS peaks (relative to basepeak) (default: %(default)s)", default=0.1,type=float)
        sc.add_argument('-i', '--ionisation', help="Ionisation mode (default: %(default)s)", default="1", choices=["-1", "1"])
        # sygma arguments
        sc.add_argument('reactants', type=argparse.FileType('r'), help="File with smiles used as reactants")
        sc.add_argument('-n', '--n_reaction_steps', help="Maximum number of reaction steps (default: %(default)s)", default=2,type=int)
        sc.add_argument('-m', '--metabolism_types', help="1 and/or 2 for phase 1 and 2 biotransformations (default: %(default)s)", nargs='+', default=["phase1", "phase2"], choices=["phase1", "phase2"])
        # output arguments
        sc.add_argument('db', type=argparse.FileType('w'), help="Sqlite database file with results")
        sc.add_argument('-b', '--max_broken_bonds', help="Maximum number of bonds broken in substructures generated from metabolites (default: %(default)s)", default=4,type=int)
        sc.add_argument('-a', '--abs_peak_cutoff', help="abs intensity threshold for storing peaks in database (default: %(default)s)", default=1000,type=float)
        sc.add_argument('-r', '--rel_peak_cutoff', help="fraction of basepeak intensity threshold for storing peaks in database (default: %(default)s)", default=0.01,type=float)
        sc.add_argument('--precursor_mz_precision', help="precision for matching precursor mz with peak mz in parent scan (default: %(default)s)", default=0.001,type=float)
        sc.add_argument('-u', '--use_msms_only', help="annotate only peaks with fragmentation data (default: %(default)s)", default=True)
        sc.add_argument('-f', '--use_fragmentation', default=True)
 
        sc.set_defaults(func=self.allinone)

        sc = subparsers.add_parser("metabolize", help=self.metabolize.__doc__, description=self.metabolize.__doc__)
        sc.add_argument('db', type=argparse.FileType('r'), help="Sqlite database file with results")
        # sygma arguments
        sc.add_argument('reactants', type=argparse.FileType('r'), help="File with smiles used as reactants")
        sc.add_argument('-n', '--n_reaction_steps', help="Maximum number of reaction steps (default: %(default)s)", default=2)
        sc.add_argument('-m', '--metabolism_types', help="1 and/or 2 for phase 1 and 2 biotransformations (default: %(default)s)", default="12", choices=["1", "2", "12"])
        # output arguments
        sc.set_defaults(func=self.metabolize)

        sc = subparsers.add_parser("mzxml", help=self.mzxml.__doc__, description=self.mzxml.__doc__)
        # mzxml arguments
        sc.add_argument('mzxml', type=argparse.FileType('r'), help="mzXMl file with MS/MS data")
        sc.add_argument('-p', '--mz_precision', help="precision in Dalton (default: %(default)s)", default=0.01)
        sc.add_argument('-c', '--msfilter', help="cutoff value to filter MS peaks (absolute) (default: %(default)s)", default=2e5)
        sc.add_argument('-d', '--msmsfilter', help="cutoff value to filter MSMS peaks (relative to basepeak) (default: %(default)s)", default=0.1)
        sc.add_argument('-i', '--ionisation', help="Ionisation mode (default: %(default)s)", default="1", choices=["-1", "1"])
        # output arguments
        sc.add_argument('db', type=argparse.FileType('w'), help="Sqlite database file with results")
        sc.set_defaults(func=self.mzxml)

        sc = subparsers.add_parser("sd2smiles", help=self.sd2smiles.__doc__, description=self.sd2smiles.__doc__)
        sc.add_argument('input', type=argparse.FileType('r'), help="Sd file")
        sc.add_argument('output', type=argparse.FileType('w'), help="File with smiles which can be used as metabolite reactantss")
        sc.set_defaults(func=self.sd2smiles)

        sc = subparsers.add_parser("smiles2sd", help=self.smiles2sd.__doc__, description=self.smiles2sd.__doc__)
        sc.add_argument('input', type=argparse.FileType('r'), help="File with smiles which can be used as metabolite reactantss")
        sc.add_argument('output', type=argparse.FileType('w'), help="Sd file")
        sc.add_argument('-d', '--depiction', help="Depiction (default: %(default)s)", default="1D", choices=["1D", "2D", "3D"])
        sc.set_defaults(func=self.smiles2sd)

    def version(self):
        return '1.0' # TODO move to main magma package and reuse in setup.py so version is specified in one place

    def allinone(self, args):
        """Reads reactants file and MS/MS datafile, generates metabolites from reactants and matches them to peaks"""
        magma.set_DB(args.db.name)
        for mol in self.smiles2mols(args.reactants):
            magma.add_metabolite(
                                 Chem.MolToMolBlock(mol).encode('ascii'),
                                 mol.GetProp('_Name'),
                                 1.0, 0, 'PARENT', 1)
        magma.metabolize_all(args.metabolism_types, args.n_reaction_steps)
        magma.commit_DB()
        magma.set_run_data(
                           n_reaction_steps=args.n_reaction_steps,
                           metabolism_types=args.metabolism_types,
                           ms_filename=args.mzxml.name,
                           ionisation=args.ionisation,
                           use_fragmentation=args.use_fragmentation,
                           max_broken_bonds=args.max_broken_bonds,
                           abs_peak_cutoff=args.abs_peak_cutoff,
                           rel_peak_cutoff=args.rel_peak_cutoff,
                           mz_precision=args.mz_precision,
                           precursor_mz_precision=args.precursor_mz_precision,
                           ms_intensity_cutoff=args.ms_intensity_cutoff,
                           msms_intensity_cutoff=args.msms_intensity_cutoff,
                           use_msms_only=args.use_msms_only
                           )        
        magma.storeMZxmlFile(args.mzxml.name)
        magma.commit_DB()
        magma.buildspectra()
        magma.searchAllMetabolites()
        magma.commit_DB()

    def metabolize(self, args):
        """Reads reactants file and existing result database, generates metabolites from reactants and matches them to peaks"""
        magma.set_DB(args.db.name)
        for mol in self.smiles2mols(args.reactants):
            magma.add_metabolite(
                                 Chem.MolToMolBlock(mol),
                                 mol.GetProp('_Name'),
                                 1.0, 0, 'PARENT', 1)
        magma.metabolize_all(args.metabolism_types, args.nsteps)
        magma.commit_DB()

    def mzxml(self, args):
        """Reads MS/MS datafile"""
        magma.set_DB(args.db.name)
        magma.set_run_data(None, None, None, args.mzxml.name, args.ionisation, True, args.msfilter, args.msmsfilter, args.precision, True)
        magma.storeMZxmlFile(args.mzxml)
        magma.commit_DB()

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

    def smiles2mols(self, smiles):
        mols = []
        for line in smiles:
            line = line.strip()
            if line=="":
                continue
            (smilestring, molname) = line.split('|')
            mol = Chem.MolFromSmiles(smilestring)
            mol.SetProp('_Name', molname.encode('ascii','xmlcharrefreplace'))
            AllChem.Compute2DCoords(mol)
            mols.append(mol)
        return mols

    def run(self):
        """Parse arguments and runs subcommand"""
        args = self.parser.parse_args()
        return args.func(args)
