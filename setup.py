import os
import sys
import glob

from setuptools import setup, find_packages

setup(
	name = "perovgen",
	version = "3.6.4",
	license = "Hanbat National Univ.",
	description = "Automated VASP Program with Pymatgen",
	author = "Jong Goo Park",
	author_email = "jgp505@gmail.com",
	#url = "",
	#downloadurl = "",
	install_requires= ["numpy>=1.9","pymatgen>=2022.0.7","matplotlib>=1.1","palettable>=3.2.0"],
	packages = find_packages(),
	python_requires = ">=3",
	zip_safe = False,
	clssifiers = [
		'Programming Language :: Python :: 3',
		'Programming Language :: Python :: 3.2',
		'Programming Language :: Python :: 3.3',
		'Programming Language :: Python :: 3.4',
		'Programming Language :: Python :: 3.5',
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: 3.7',\
		],
	scripts = glob.glob(os.path.join(os.path.dirname(__file__),"argument",'*'))
	)

