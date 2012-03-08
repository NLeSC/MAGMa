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

        sc = subparsers.add_parser("all_in_one", help=self.all_in_one.__doc__, description=self.all_in_one.__doc__)
        sc.add_argument('-z', '--description', help="Description of the job (default: %(default)s)", default="",type=str)
        # read_ms_data arguments
        sc.add_argument('ms_data', type=argparse.FileType('r'), help="file with MS/MS data")
        sc.add_argument('--ms_data_format', help="MS data input format (default: %(default)s)", default="mzxml", choices=["mzxml"])
        sc.add_argument('-l', '--max_ms_level', help="Maximum ms level to be processsed (default: %(default)s)", default=10,type=int)
        sc.add_argument('-a', '--abs_peak_cutoff', help="abs intensity threshold for storing peaks in database (default: %(default)s)", default=1000,type=float)
        sc.add_argument('-r', '--rel_peak_cutoff', help="fraction of basepeak intensity threshold for storing peaks in database (default: %(default)s)", default=0.01,type=float)
        # add_structures arguments
        sc.add_argument('structures', type=argparse.FileType('rb'), help="File with smiles used as structures")
        sc.add_argument('-t', '--structure_format', help="Structure input type (default: %(default)s)", default="smiles", choices=["smiles", "sdf"])
        sc.add_argument('-s', '--n_reaction_steps', help="Maximum number of reaction steps (default: %(default)s)", default=2,type=int)
        sc.add_argument('-m', '--metabolism_types', help="1 and/or 2 for phase 1 and 2 biotransformations (default: %(default)s)", default="phase1,phase2", type=str)
        # annotate arguments
        sc.add_argument('-p', '--mz_precision', help="precision in Dalton (default: %(default)s)", default=0.001,type=float)
        sc.add_argument('-c', '--ms_intensity_cutoff', help="cutoff value to filter MS peaks (absolute) (default: %(default)s)", default=1e6,type=float)
        sc.add_argument('-d', '--msms_intensity_cutoff', help="cutoff value to filter MSMS peaks (relative to basepeak) (default: %(default)s)", default=0.1,type=float)
        sc.add_argument('-i', '--ionisation_mode', help="Ionisation mode (default: %(default)s)", default="1", choices=["-1", "1"])
        sc.add_argument('-b', '--max_broken_bonds', help="Maximum number of bonds broken in substructures generated from metabolites (default: %(default)s)", default=4,type=int)
        sc.add_argument('--precursor_mz_precision', help="precision for matching precursor mz with peak mz in parent scan (default: %(default)s)", default=0.005,type=float)
        sc.add_argument('-u', '--use_all_peaks', help="annotate all peaks, including those without fragmentation data (default: %(default)s)", action="store_true")
        sc.add_argument('-f', '--skip_fragmentation', action="store_true")
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.all_in_one)

        sc = subparsers.add_parser("add_structures", help=self.add_structures.__doc__, description=self.add_structures.__doc__)
        sc.add_argument('-z', '--description', help="Description of the job (default: %(default)s)", default="",type=str)
        # add_structures arguments
        sc.add_argument('structures', type=argparse.FileType('rb'), help="File with smiles used as structures")
        sc.add_argument('-t', '--structure_format', help="Structure input type (smiles or sdff) (default: %(default)s)", default="smiles", choices=["smiles", "sdf"])
        sc.add_argument('-s', '--n_reaction_steps', help="Maximum number of reaction steps (default: %(default)s)", default=2,type=int)
        sc.add_argument('-m', '--metabolism_types', help="1 and/or 2 for phase 1 and 2 biotransformations (default: %(default)s)", default="phase1,phase2", type=str)
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.add_structures)

        sc = subparsers.add_parser("read_ms_data", help=self.read_ms_data.__doc__, description=self.read_ms_data.__doc__)
        sc.add_argument('-z', '--description', help="Description of the job (default: %(default)s)", default="",type=str)
        # read_ms_data arguments
        sc.add_argument('ms_data', type=argparse.FileType('r'), help="file with MS/MS data")
        sc.add_argument('--ms_data_format', help="MS data input format (default: %(default)s)", default="mzxml", choices=["mzxml"])
        sc.add_argument('-l', '--max_ms_level', help="Maximum ms level to be processsed (default: %(default)s)", default=10,type=int)
        sc.add_argument('-a', '--abs_peak_cutoff', help="abs intensity threshold for storing peaks in database (default: %(default)s)", default=1000,type=float)
        sc.add_argument('-r', '--rel_peak_cutoff', help="fraction of basepeak intensity threshold for storing peaks in database (default: %(default)s)", default=0.01,type=float)
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.read_ms_data)

        sc = subparsers.add_parser("annotate", help=self.annotate.__doc__, description=self.annotate.__doc__)
        sc.add_argument('-z', '--description', help="Description of the job (default: %(default)s)", default="",type=str)
        # annotate arguments
        sc.add_argument('-p', '--mz_precision', help="precision in Dalton (default: %(default)s)", default=0.001,type=float)
        sc.add_argument('-c', '--ms_intensity_cutoff', help="cutoff value to filter MS peaks (absolute) (default: %(default)s)", default=1e6,type=float)
        sc.add_argument('-d', '--msms_intensity_cutoff', help="cutoff value to filter MSMS peaks (relative to basepeak) (default: %(default)s)", default=0.1,type=float)
        sc.add_argument('-i', '--ionisation_mode', help="Ionisation mode (default: %(default)s)", default="1", choices=["-1", "1"])
        sc.add_argument('-b', '--max_broken_bonds', help="Maximum number of bonds broken in substructures generated from metabolites (default: %(default)s)", default=4,type=int)
        sc.add_argument('--precursor_mz_precision', help="precision for matching precursor mz with peak mz in parent scan (default: %(default)s)", default=0.005,type=float)
        sc.add_argument('-u', '--use_all_peaks', help="annotate all peaks, including those without fragmentation data (default: %(default)s)", action="store_true")
        sc.add_argument('-f', '--skip_fragmentation', action="store_true")
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.annotate)

        sc = subparsers.add_parser("sd2smiles", help=self.sd2smiles.__doc__, description=self.sd2smiles.__doc__)
        sc.add_argument('input', type=argparse.FileType('rb'), help="Sd file")
        sc.add_argument('output', type=argparse.FileType('w'), help="File with smiles which can be used as metabolite reactantss")
        sc.set_defaults(func=self.sd2smiles)

        sc = subparsers.add_parser("smiles2sd", help=self.smiles2sd.__doc__, description=self.smiles2sd.__doc__)
        sc.add_argument('input', type=argparse.FileType('r'), help="File with smiles which can be used as metabolite reactantss")
        sc.add_argument('output', type=argparse.FileType('w'), help="Sd file")
        sc.add_argument('-d', '--depiction', help="Depiction (default: %(default)s)", default="2D", choices=["1D", "2D", "3D"])
        sc.set_defaults(func=self.smiles2sd)

    def version(self):
        return '1.0' # TODO move to main magma package and reuse in setup.py so version is specified in one place

    def all_in_one(self, args):
        """Reads reactants file and MS/MS datafile, generates metabolites from reactants and matches them to peaks"""

        magma_session = magma.MagmaSession(args.db,args.description)
        self._add_structures(args, magma_session)
        self._read_ms_data(args, magma_session)
        self._annotate(args, magma_session)

    def add_structures(self, args):
        """Reads reactants file and existing result database, generates metabolites from reactants and matches them to peaks"""
        magma_session = magma.MagmaSession(args.db,args.description)
        self._add_structures(args, magma_session)

    def read_ms_data(self, args):
        magma_session = magma.MagmaSession(args.db,args.description)
        self._read_ms_data(args, magma_session)

    def annotate(self, args):
        magma_session = magma.MagmaSession(args.db,args.description)
        self._annotate(args, magma_session)

    def _add_structures(self, args, magma_session):
        struct_engine = magma_session.get_structure_engine(args.metabolism_types, args.n_reaction_steps) # TODO remove arguments
        if args.structure_format == 'smiles':
            for mol in self.smiles2mols(args.structures):
                struct_engine.add_structure(Chem.MolToMolBlock(mol), mol.GetProp('_Name'), 1.0, 0, 'PARENT', 1)
        elif args.structure_format == 'sdf':
            for mol in Chem.SDMolSupplier(args.structures.name):
                struct_engine.add_structure(Chem.MolToMolBlock(mol), mol.GetProp('_Name'), 1.0, 0, 'PARENT', 1)
        struct_engine.metabolize_all(args.metabolism_types, args.n_reaction_steps)

    def _read_ms_data(self, args, magma_session):
        ms_data_engine = magma_session.get_ms_data_engine(abs_peak_cutoff=args.abs_peak_cutoff, 
            rel_peak_cutoff=args.rel_peak_cutoff,max_ms_level=args.max_ms_level)
        if args.ms_data_format == "mzxml":
            ms_data_engine.store_mzxml_file(args.ms_data.name)
        elif args.ms_data_format == "peak_list":
            pass

    def _annotate(self, args, magma_session):
        annotate_engine = magma_session.get_annotate_engine(ionisation_mode=args.ionisation_mode, 
            skip_fragmentation=args.skip_fragmentation, 
            max_broken_bonds=args.max_broken_bonds, 
            ms_intensity_cutoff=args.ms_intensity_cutoff, 
            msms_intensity_cutoff=args.msms_intensity_cutoff, 
            mz_precision=args.mz_precision, 
            precursor_mz_precision=args.precursor_mz_precision, 
            use_all_peaks=args.use_all_peaks)
        annotate_engine.build_spectra()
        annotate_engine.search_all_structures()

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
        w = Chem.SDWriter(args.output.name)
        for line in args.input:
            line = line.strip()
            if line=="":
                continue
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
            mol.SetProp('_Name', molname)
            AllChem.Compute2DCoords(mol)
            mols.append(mol)
        return mols

    def run(self):
        """Parse arguments and runs subcommand"""
        args = self.parser.parse_args()
        return args.func(args)


if __name__ == "__main__":
    main()
