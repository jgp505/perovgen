import os
import sys

import argparse
import textwrap
import itertools

from perovgen.pygmd.cli.gmd_config import config
from perovgen.pygmd.cli.gmd_plot import plot
from perovgen.pygmd.cli.gmd_replace import randomreplace
from perovgen.pygmd.cli.gmd_auto import analyze_calculation
from perovgen.pygmd.cli.gmd_analysis import analysis

__author__ = "Green Materials Design Lab"
__maintainer__ = "Park Jong Goo"
__email__ = "jgp505@gmail.com"
__version__ = "3.6.4"

def main() :
    parser = argparse.ArgumentParser(
        description='''
	gmd3 is a convenient script that uses perovgen that uses 
	the pymatgen(https://pymatgen.org/) modules to perform 
	many analyses, plotting and format conversions. Also, 
	This script used exclusively in the Green Materials Design Lab Team. 
	We research Organic-Inorganic hybrid halide Perovskite Materials. 
    This script works based on several sub-commands with their own options. 
    To see the options for the sub-commands, type “gmd3 sub-command -h”. 
	''',
        epilog = '''Version : {}'''.format(__version__)
    )
    subparsers = parser.add_subparsers()

    # Config Option
    parser_config = subparsers.add_parser("config",
    help = "Register and manage shell scripts")
    parser_config.add_argument("-s","--shell",dest="shell",
    action = "store",
    nargs = 1,
    help = "Designate the shell script required Perovgen autocalculation.")
    parser_config.add_argument("-c","--check",dest="check",
    action = "store_true",
    help = "Check the Shell Script")
    parser_config.add_argument("-r","--remove",dest='remove',
    action = 'store_true',
    help="Remove the Shell Script registered")
    parser_config.add_argument("-g","--generate",dest='generate',
    action = 'store',
    nargs = 1,
    default = False,
    help="Generated shell")
    parser_config.set_defaults(func=config)

    # Plotting Option
    parser_plot = subparsers.add_parser("plot",
    help = ''' Plotting for electronic structure Analysis through VASP program
    ''')

    groups = parser_plot.add_mutually_exclusive_group(required=True)
    groups.add_argument("--band",dest="band",
	action="store_true",
	help="Plotting band structure")
    groups.add_argument("--dos",dest="dos",
	action="store_true",
	help="Plotting DOS")
    groups.add_argument("--bdos",dest="bdos",
	action="store_true",
    help="Plotting both Band structure and DOS")
    groups.add_argument("--pdos",dest="partial",
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
    Substitutes a composition from element to element or element to molecule.

    1. Element => Element
    2. Element => Molecule
    3. Relaxation structure (beta test)
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

    #parser_mole.add_argument("-i","--input",dest="input",type=str,
    #default=False,
    #help="read the input file(input.gmd). Input is composed of STRUC, INCAR, METHOD, SHELL, KPOINTS")
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
    parser_auto.add_argument("-i",dest="inputpath",
    help="read the input.gmd file. Input is composed of INCAR, CALMODE, SHELL, KPOINTS")
    parser_auto.add_argument('-p',"--path",dest="strucpath",
    nargs="*",action = "store",help=" \n\nSpecifies the path of the structural file or files Default is to bring the STRUC in the inputfile(input.gmd)\n")

    parser_auto.add_argument("-d","--directory",dest="directory",
    action="store_true",
    default=False,
    help="An optional argument to enter when you don't want to run the script")
    parser_auto.add_argument("--de",dest="deleteselect",action="store_true",
    default=False, help="remove the selective dynamic from the POSCAR")
    parser_auto.add_argument("--soc",dest="soc",action="store_true",default=False,
    help='''
    Calculate Spin Orbit Coupling (SOC)
    ''')
    parser_auto.set_defaults(func=analyze_calculation)

    ### Analysis
    parser_analysis = subparsers.add_parser("analysis",
    help = '''
    Support analysis  
    ''')
    parser_analysis.add_argument("-g",'--graph.yaml',dest='graph',
    choices = ['BD','B','D'],
    help = "Generate the graph.yaml that plot the bands structure.")
    parser_analysis.add_argument("-i","--input",dest='input',
    action='store_true',
    default=False,
    help='generated input.gmd')
    parser_analysis.add_argument("-pdos",dest='pdos',
    action='store_true',
    default=False,
    help='do pdos')
    parser_analysis.add_argument("-c","--convert",dest='convert',
    action='store',
    nargs='*',
    default=False,
    help='convert POSCAR to cif')
    parser_analysis.set_defaults(func=analysis)

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
