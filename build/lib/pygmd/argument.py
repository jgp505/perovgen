import os
import sys

import argparse
import textwrap
import itertools

from perovgen.pygmd.cli.gmd_config import config
from perovgen.pygmd.cli.gmd_files import files
from perovgen.pygmd.cli.gmd_plot import plot
from perovgen.pygmd.cli.gmd_replace import randomreplace
from perovgen.pygmd.cli.gmd_auto import analyze_calculation
from perovgen.pygmd.cli.gmd_data import dataframe

__author__ = "Green Materials Design Lab"
__maintainer__ = "Park Jong Goo"
__email__ = "jgp505@gmail.com"
__version__ = "3.5.1"

def main() :
    parser = argparse.ArgumentParser(
        description='''
	gmd is a convenient script that uses gmdkit that uses 
	the pymatgen(https://pymatgen.org/) modules to perform 
	many analyses, plotting and format conversions. Also, 
	This script used exclusively in the GMD Lab. 
	Because we research Organic-Inorganic hybrid halide Perovskite Materials, 
	we focus on the that materials. This script works based on several 
	sub-commands with their own options. To see the options for the sub-commands, type “gmd  sub-command -h”. 
	''',
        epilog = '''Version : {}'''.format(__version__)
    )
    subparsers = parser.add_subparsers()

    # Config Option
    parser_config = subparsers.add_parser("config",
    help = "Manage shell script for VASP calculation")
    parser_config.add_argument("-s","--shell",dest="shell",
    action = "store",
    nargs = 1,
    help = "Designate the shell script required the VASP calculation.")
    parser_config.add_argument("-c","--check",dest="check",
    action = "store_true",
    help = "Check the Shell Script")
    parser_config.add_argument("-r","--remove",dest='remove',
    action = 'store_true',
    help="Remove the Shell Script")
    parser_config.set_defaults(func=config)

    # File Option
    parser_file = subparsers.add_parser("file",
    help = ''' Structure file for VASP calculation (e.g. POSCAR, CONTCAR)
    and file for electronic structure plotting ''')
    parser_file.add_argument("-g","--graph", dest="graph", 
    choices=["BD","B","D","All"],
    help = "Generate the graph.yaml that plot the bands structure.")	
    parser_file.add_argument("-s","--strain",dest="strain", 
    action = "store", 
    nargs=2,
    help="Deformate from structure file to make the modified structure file.")
    parser_file.add_argument("-i","--input",dest="input",
    action = "store",
    nargs="*",
    default = False,
    help="The structure files in the current directory are printed")
    parser_file.add_argument("-pdos",dest='pdos',
    action='store_true',
    help='make dos file with pdos')
    parser_file.add_argument("--mp",dest='mp',
    action='store',
    nargs='*',
    help='TEST')
    parser_file.set_defaults(func=files)

    # Plotting Option
    parser_plot = subparsers.add_parser("plot",
    help = ''' Plotting for electronic structure Analysis through VASP program
    (e.g. Band Structure, DOS)''')

    groups = parser_plot.add_mutually_exclusive_group(required=True)
    groups.add_argument("-b","--band",dest="band",
	action="store_true",
	help="Plotting band structure")
    groups.add_argument("-d","--dos",dest="dos",
	action="store_true",
	help="Plotting DOS")
    groups.add_argument("-bd","--bdos",dest="bdos",
	action="store_true",
    help="Plotting both Band structure and DOS")
    groups.add_argument("-p","--partial",dest="partial",
    action="store_true",
    help="Plotting p-DOS that required LIST file")
    groups.add_argument("-bc","--check",dest="bandcheck",
	action = "store_true",
    help = "Check the Band Gap, CBM and VBM points")
    
    parser_plot.add_argument("-n","--name", dest="name",
    action="store",
    nargs="*",
    default="dos",
    help="Bring the dos file made with pdos")
    parser_plot.add_argument("--path",dest="path",
    default="vasprun.xml", nargs=1,
    help="Bring the vasprun.xml.default is to bring the vasprun.xml in current folder")
    parser_plot.add_argument("--format",dest="format",
    default="pdf", nargs=1,
    help="It changes the format of the file to be saved with --outfile")
    parser_plot.add_argument("--outfile", dest="out_file", 
    type=str,
	help="Save plot to file instead of displaying.")
    parser_plot.set_defaults(func=plot)

    # Molecule 
    parser_mole = subparsers.add_parser("replace",
    help = '''
    In the Inorganic Perovksite(e.g. CsPbBr3), the inorganic material
    is substituted with an organic material or elements
    (e.g. MA(Methylammonium),FA(Formamidinium)) 
    ''')
    groups = parser_mole.add_mutually_exclusive_group(required=True)
    groups.add_argument("--ma",dest="ma",
    action = "store_true",
    help = "Substitute with MA")
    groups.add_argument("--fa",dest="fa",
    action = "store_true",
    help = "Substitute with FA")
    groups.add_argument("--gua",dest="gua",
    action = "store_true",
    help = "Substitute with GUA")
    groups.add_argument("--dima",dest="dima",
    action = "store_true",
    help = "Substitute with dimethylammonium")
    groups.add_argument("--trima",dest="trima",
    action = "store_true",
    help = "Substitute with trimethylammonium")
    groups.add_argument("--tetrama",dest="tetrama",
    action = "store_true",
    help = "Substitute with tetramethylammonium")
    groups.add_argument("--zolium",dest="zolium",
    action = "store_true",
    help = "Substitute with Imidazolium")
    groups.add_argument("--ele",dest="sub",
    action = "store_true",
    help="Substitute with Element")

    parser_mole.add_argument("--path", dest="path",
    type=str,
    default = os.getcwd(),
    nargs="*",
    action = "store",
    help = "The structural of file format or this folder in files will be read")
    parser_mole.add_argument("--degree","-d",dest="degree",
    action = "store_false",default = True,
    help = "Specify the Angle of Tilting of MA and FA")
    parser_mole.add_argument("--position","-p",dest="position",
    action = "store_false",default = True,
    help = "Specify the Position of Structure that loaded name")
    parser_mole.add_argument("--csv","-c",dest="csv",
    action = "store_false", default= True,
    help = "Save the molecule information for CSV file")

    parser_mole.add_argument("-i","--input",dest="input",type=str,
    default=False,
    help="read the input file(input.gmd). Input is composed of STRUC, INCAR, METHOD, SHELL, KPOINTS")
    parser_mole.set_defaults(func=randomreplace)

    # autocalculation
    parser_auto = subparsers.add_parser("auto",
    help = '''
    As argument for vasp automatic calculation
    The kind of calculation is [PBE,PBESOL,VDW,SCAN] and Default
    Read the structure file and input file necessary
    for vasp calculation automatically.
    If shell script is designated through config, 
    vasp calculation can be proceed with the [--run] command
    ''')
    parser_auto.add_argument("-i","--input",dest="input",type=str,
    help="read the input file(input.gmd). Input is composed of STRUC, INCAR, METHOD, SHELL, KPOINTS")
    parser_auto.add_argument("-p","--path",dest="path",
    type=str,
    default = False,
    nargs="*",
    action = "store",
    help='''
    \n\nSpecifies the path of the structural file or files
    Default is to bring the STRUC in the inputfile(input.gmd)\n
    ''')
    parser_auto.add_argument("-d","--directory",dest="directory",
    action="store_true",
    default=False,
    help="An optional argument to enter when you don't want to run the script")
    parser_auto.add_argument("--de",dest="deleteselect",action="store_true",
    default=False, help="remove the selective dynamic from the POSCAR")
    parser_auto.add_argument("--soc",dest="soc",action="store_true",default=False,
    help='''
    If you want to calculate SOC(Spin Orbit Coupling), 
    Use -o or --orbit mode insted of -s or --structure
    If you want to run calculation, adding --run [SHELL SCRIPT NAME] [NAMINNG] at the end.
    But, NAMING can be omitted.''')
    parser_auto.set_defaults(func=analyze_calculation)

    # Dataframe
    parser_data = subparsers.add_parser("data",
    help = '''
    By designating the folder where the calculation is completed,
    and the structured data is saved as a csv file 
    according to the specified platform.
    ''')
    parser_data.add_argument("-p","--path",dest="path",
    type=str,
    default = os.getcwd(),
    nargs="*",
    action = "store",
    help='''
    \n\nDesignated path\n\n
    ''')
    parser_data.add_argument("-b","--band",dest="band",
    action="store_true",
    help="Add the fermi energy and band gap")
    parser_data.add_argument("-n","--name",dest="name",
    default = "Default",
    type = str,
    help="csv file name")
    parser_data.set_defaults(func=dataframe)

    try :
        import argcomplete
        argcomplete.autocomplete(parser)
    except ImportError :
        pass
    args = parser.parse_args()
    try :
        getattr(args,"func")
    except AttributeError :
        parser.print_help()
        sys.exit(0)
    args.func(args)

if __name__ == "__main__" :
    main()
