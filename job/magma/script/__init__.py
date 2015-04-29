import argparse
import sys,logging,shutil
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

        sc = subparsers.add_parser("init_db", help=self.init_db.__doc__, description=self.init_db.__doc__)
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.init_db)

        sc = subparsers.add_parser("add_structures", help=self.add_structures.__doc__, description=self.add_structures.__doc__)
        sc.add_argument('-z', '--description', help="Description of the job (default: %(default)s)", default="",type=str)
        # add_structures arguments
        sc.add_argument('-t', '--structure_format', help="Structure input type (default: %(default)s)", default="smiles", choices=["smiles", "sdf"])
        sc.add_argument('-p', '--pubchem_names', help="Get references to PubChem (default: %(default)s)", action="store_true")
        sc.add_argument('-m', '--mass_filter', help="Filter input structures on maximum monoisotopic mass (default: %(default)s)", default=9999,type=int)
        sc.add_argument('-l', '--log', help="Set logging level (default: %(default)s)", default='info',choices=['debug','info','warn','error'])
        sc.add_argument('structures', type=str, help="File with structures, or a single smiles string")
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.add_structures)

        sc = subparsers.add_parser("metabolize", help=self.metabolize.__doc__, description=self.metabolize.__doc__)
        sc.add_argument('-z', '--description', help="Description of the job (default: %(default)s)", default="",type=str)
        # metabolize arguments
        sc.add_argument('-j', '--molids', type=argparse.FileType('rb'), help="File with structure ids")
        sc.add_argument('-n', '--n_reaction_steps', help="Maximum number of reaction steps (default: %(default)s)", default=1,type=int)
        sc.add_argument('-m', '--metabolism_types', help="digest,gut,phase1,phase2, (default: %(default)s)", default="phase1,phase2", type=str)
        sc.add_argument('-s', '--scenario', default=None, type=str, help="""Scenario file, each line defines a separate stage:
                                        action(glycosidase/gut/phase1[_selected]/phase2[_selected]/mass_filter),value(nsteps/mass limit)""")
        sc.add_argument('-t', '--time_limit', help="Maximum allowed time in minutes (default: %(default)s)", default=None,type=float)
        sc.add_argument('-p', '--pubchem_names', help="Get references to PubChem (default: %(default)s)", action="store_true")
        sc.add_argument('-l', '--log', help="Set logging level (default: %(default)s)", default='info',choices=['debug','info','warn','error'])
        sc.add_argument('--call_back_url', help="Call back url (default: %(default)s)", default=None,type=str)
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.metabolize)

        sc = subparsers.add_parser("read_ms_data", help=self.read_ms_data.__doc__, description=self.read_ms_data.__doc__)
        sc.add_argument('-z', '--description', help="Description of the job (default: %(default)s)", default="",type=str)
        # read_ms_data arguments
        sc.add_argument('ms_data', type=str, help="file with MS/MS data")
        sc.add_argument('-f', '--ms_data_format', help="MS data input format (default: %(default)s)", default="mzxml", choices=["mzxml", "mass_tree","form_tree_pos","form_tree_neg"])
        sc.add_argument('-i', '--ionisation_mode', help="Ionisation mode (default: %(default)s)", default="1", choices=["-1", "1"])
        sc.add_argument('-m', '--max_ms_level', help="Maximum MS level to be processsed (default: %(default)s)", default=10,type=int)
        sc.add_argument('-a', '--abs_peak_cutoff', help="Absolute intensity threshold for storing peaks in database (default: %(default)s)", default=1000,type=float)
        sc.add_argument('-p', '--mz_precision', help="Maximum relative m/z error (ppm) (default: %(default)s)", default=5,type=float)
        sc.add_argument('-q', '--mz_precision_abs', help="Maximum absolute m/z error (Da) (default: %(default)s)", default=0.001,type=float)
        sc.add_argument('--precursor_mz_precision', help="Maximum absolute error of precursor m/z values (default: %(default)s)", default=0.005,type=float)
        sc.add_argument('-s', '--scan', help="Read only spectral tree specified by MS1 scan number (default: %(default)s)", default=None,type=str)
        sc.add_argument('-t', '--time_limit', help="Maximum allowed time in minutes (default: %(default)s)", default=None,type=float)
        sc.add_argument('-l', '--log', help="Set logging level (default: %(default)s)", default='info',choices=['debug','info','warn','error'])
        sc.add_argument('--call_back_url', help="Call back url (default: %(default)s)", default=None,type=str)
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.read_ms_data)

        sc = subparsers.add_parser("annotate", help=self.annotate.__doc__, description=self.annotate.__doc__)
        sc.add_argument('-z', '--description', help="Description of the job (default: %(default)s)", default="",type=str)
        # annotate arguments
        sc.add_argument('-c', '--ms_intensity_cutoff', help="Minimum intensity of MS1 precursor ion peaks to be annotated (default: %(default)s)", default=1e6,type=float)
        sc.add_argument('-d', '--msms_intensity_cutoff', help="Minimum intensity of of fragment peaks to be annotated, as percentage of basepeak (default: %(default)s)", default=5,type=float)
        sc.add_argument('-j', '--molids', help="structure ids, comma separated", type=str)
        sc.add_argument('-b', '--max_broken_bonds', help="Maximum number of bond breaks to generate substructures (default: %(default)s)", default=3,type=int)
        sc.add_argument('-w', '--max_water_losses', help="Maximum number of additional water (OH) and/or ammonia (NH2) losses (default: %(default)s)", default=1,type=int)
        sc.add_argument('-u', '--use_all_peaks', help="Annotate all level 1 peaks, including those not fragmented (default: %(default)s)", action="store_true")
        sc.add_argument('--skip_fragmentation', help="Skip substructure annotation of fragment peaks (default: %(default)s)", action="store_true")
        sc.add_argument('-f', '--fast', help="Quick calculations for molecules up to 64 atoms (default: %(default)s)", action="store_true")
        sc.add_argument('-s', '--structure_database', help="Retrieve molecules from structure database  (default: %(default)s)", default="", choices=["pubchem","kegg","hmdb","metacyc"])
        sc.add_argument('-o', '--db_options', help="Specify structure database option: db_filename,max_mim,max_64atoms,incl_halo,min_refscore(only for PubChem) (default: %(default)s)",default=",1200,False",type=str)
        sc.add_argument('-a', '--adducts' ,default=None,type=str, help="""Specify adduct (as comma separated list) for matching at MS1.
                                                                        Positive mode: [Na,K,NH4] Negative mode: [Cl]
                                                                        (default: %(default)s)""")
        sc.add_argument('-m', '--max_charge', help="Maximum charge state (default: %(default)s)", default=1,type=int)
        sc.add_argument('-n', '--ncpus', help="Number of parallel cpus to use for annotation (default: %(default)s)", default=1,type=int)
        sc.add_argument('--scans', help="Search in specified scans (default: %(default)s)", default="all",type=str)
        sc.add_argument('-t', '--time_limit', help="Maximum allowed time in minutes (default: %(default)s)", default=None,type=float)
        sc.add_argument('-l', '--log', help="Set logging level (default: %(default)s)", default='info',choices=['debug','info','warn','error'])
        sc.add_argument('--call_back_url', help="Call back url (default: %(default)s)", default=None,type=str)
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.annotate)

        sc = subparsers.add_parser("select", help=self.select.__doc__, description=self.select.__doc__)
        sc.add_argument('-f', '--frag_id', help="Fragment_identifier as selection query (default: %(default)s)", default=None,type=int)
        sc.add_argument('db_in', type=str, help="Input sqlite database file with annotation results")
        sc.add_argument('db_out', type=str, help="Output qlite database file with selected results")
        sc.set_defaults(func=self.select)

        sc = subparsers.add_parser("export_structures", help=self.export_structures.__doc__, description=self.export_structures.__doc__)
        sc.add_argument('-f', '--filename', help="Output filename (default: stdout)", default=None, type=str)
        sc.add_argument('-a', '--assigned', help="Only assigned molecules (default: %(default)s)", action="store_true")
        sc.add_argument('db', type=str, help="Sqlite database file with results")
        sc.set_defaults(func=self.export_structures)
        
    def version(self):
        return '1.0' # TODO move to main magma package and reuse in setup.py so version is specified in one place

    def get_magma_session(self, db, description="",log='warn'):
        return magma.MagmaSession(db, description,log)

    def init_db(self,args):
        """Initialize database"""
        magma_session = self.get_magma_session(args.db,"")

    def add_structures(self, args, magma_session=None):
        try:
            if magma_session == None:
                magma_session = self.get_magma_session(args.db,args.description,args.log)
            struct_engine = magma_session.get_structure_engine(pubchem_names=args.pubchem_names)
            if args.structure_format == 'smiles':
                struct_engine.read_smiles(args.structures,mass_filter=args.mass_filter)
            elif args.structure_format == 'sdf':
                struct_engine.read_sdf(args.structures,args.mass_filter)
            magma_session.commit()
        except Exception as error:
            if args.log == 'debug':
                logging.exception(error) #This includes a traceback
            else:
                logging.error(error) #This only prints an error message

    def metabolize(self, args, magma_session=None):
        try:
            if magma_session == None:
                magma_session = self.get_magma_session(args.db,args.description,args.log)
            struct_engine = magma_session.get_structure_engine(pubchem_names=args.pubchem_names, call_back_url=args.call_back_url)
            if args.scenario != None:
                # Convert comma separated file into list of [action,value] pairs: [[action,value][action,value]..]
                scenario=[]
                scenario_file=open(args.scenario,'r')
                for line in scenario_file:
                    step=line.split('#')[0].rstrip().split(',') # Comments indicated by # are allowed
                    if len(step)>1:
                        scenario.append(step)
                struct_engine.run_scenario(scenario,args.time_limit)
            elif args.molids == None:
                molids=struct_engine.metabolize_all(args.metabolism_types, args.n_reaction_steps)
                for molid in molids:
                    print molid
            else:
                molids=[]
                for reactantid in args.molids:
                    molids.extend(struct_engine.metabolize(reactantid,args.metabolism_types, args.n_reaction_steps))
                for molid in set(molids):
                    print molid
            magma_session.fill_molecules_reactions()
        except Exception as error:
            if args.log == 'debug':
                logging.exception(error)
            else:
                logging.error(error)

    def read_ms_data(self, args, magma_session=None):
        try:
            if magma_session == None:
                magma_session = self.get_magma_session(args.db,args.description,args.log)
            ms_data_engine = magma_session.get_ms_data_engine(ionisation_mode=args.ionisation_mode,
                    abs_peak_cutoff=args.abs_peak_cutoff,
                    mz_precision=args.mz_precision,
                    mz_precision_abs=args.mz_precision_abs,
                    precursor_mz_precision=args.precursor_mz_precision,
                    max_ms_level=args.max_ms_level,
                    call_back_url=args.call_back_url)
            if args.ms_data_format == "mzxml":
                ms_data_engine.store_mzxml_file(args.ms_data,args.scan,args.time_limit)
            else:
                tree_type={"mass_tree":0,"form_tree_neg":-1,"form_tree_pos":1}[args.ms_data_format]
                ms_data_engine.store_manual_tree(args.ms_data,tree_type)
        except Exception as error:
            if args.log == 'debug':
                logging.exception(error)
            else:
                logging.error(error)

    def annotate(self, args, magma_session=None):
        try:
            if magma_session == None:
                magma_session = self.get_magma_session(args.db,args.description,args.log)
            annotate_engine = magma_session.get_annotate_engine(skip_fragmentation=args.skip_fragmentation,
                max_broken_bonds=args.max_broken_bonds,
                max_water_losses=args.max_water_losses,
                ms_intensity_cutoff=args.ms_intensity_cutoff,
                msms_intensity_cutoff=args.msms_intensity_cutoff,
                use_all_peaks=args.use_all_peaks,
                adducts=args.adducts,
                max_charge=args.max_charge,
                call_back_url=args.call_back_url)
            if args.scans == 'all':
                scans='all'
            else:
                scans=set([])
                for s in args.scans.split(','):
                   scans.add(int(s))
            annotate_engine.build_spectra(scans)
            pubchem_molids=[]
            if args.structure_database != "":
                db_opts=['','','','','']
                db_options=args.db_options.split(',')
                for x in range(len(db_options)):
                    db_opts[x]=db_options[x]
                if args.structure_database == 'pubchem':
                    query_engine=magma.PubChemEngine(db_opts[0],(db_opts[2]=='True'),db_opts[3],db_opts[4])
                elif args.structure_database == 'kegg':
                    query_engine=magma.KeggEngine(db_opts[0],(db_opts[2]=='True'),db_opts[3])
                elif args.structure_database == 'hmdb':
                    query_engine=magma.HmdbEngine(db_opts[0],(db_opts[2]=='True'))
                elif args.structure_database == 'metacyc':
                    query_engine=magma.MetaCycEngine(db_opts[0],(db_opts[2]=='True'))
                pubchem_molids=annotate_engine.get_db_candidates(query_engine,db_opts[1])
            if args.molids == None:
                annotate_engine.search_structures(ncpus=args.ncpus,fast=args.fast,time_limit=args.time_limit)
            else:
                molids=args.molids.split(',')+pubchem_molids
                annotate_engine.search_structures(molids=molids,ncpus=args.ncpus,fast=args.fast,time_limit=args.time_limit)
            magma_session.commit()
            magma_session.fill_molecules_reactions()
                # annotate_engine.search_some_structures(molids)
        except Exception as error:
            if args.log == 'debug':
                logging.exception(error)
            else:
                logging.error(error)

    def select(self, args):
        shutil.copy(args.db_in,args.db_out)
        magma_session = self.get_magma_session(args.db_out)
        select_engine = magma_session.get_select_engine()
        select_engine.select_fragment(args.frag_id)

    def export_structures(self, args, magma_session=None):
        if magma_session == None:
            magma_session = self.get_magma_session(args.db)
        export_engine = magma_session.get_export_molecules_engine()
        if args.assigned:
            export_engine.export_assigned_molecules(args.filename)
        else:
            export_engine.export_molecules(args.filename)

    def run(self, argv=sys.argv[1:]):
        """Parse arguments and run subcommand"""
        args = self.parser.parse_args(argv)
        return args.func(args)

if __name__ == "__main__":
    main()
