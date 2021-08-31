# perovgen

perovgen is a convenient script that uses gmdkit that uses
the pymatgen(https://pymatgen.org/) modules to perform 
many analyses, plotting and format conversions.  
This script works based on several sub-commands with their own options. 
To see the options for the sub-commands, type “gmd  sub-command -h”. 

# Installation

1. Install git, if not already packaged with your system.
2. Download the perovgen code using under line :

	> git clone https://github.com/jgp505/python-perovgen.git
	
3. Navigate to perovget root direcotry :

	> cd perovgen
	
4. Install the code, using the command :
5. Configure your VASP pseudopotentails for use with pymatgen ::
	> pmg config -p <EXTRACTED_VASP_POTCAR> <MY_PSP>

	> pmg config --add PMG_VASP_PSP_DIR <MY_PSP>

Here <EXTRACTED_VASP_POTCAR> is the folder where your pseudopotentials are present and <MY_PSP> is the directory where the layout of pseudopotentials is organized by pymatgen. For more information refer to pymatgen installation instructions.


# How to use

> gmd3 -h

# Authors 

Jong-Goo Park - Owner of this project

# License

This project is licensed under the MIT License - see the LICENSE file for details

# Links

>[Github](https://github.com/jgp505/python-perovgen.git)

>[pip](https://pypi.org/manage/project/perovgen/releases)
