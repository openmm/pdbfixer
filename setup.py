"""pdbfixer Fixes problems in PDB files

"""
from __future__ import print_function
import os
from os.path import relpath, join
from setuptools import setup, find_packages
DOCLINES = __doc__.split("\n")

########################
__version__ = '1.0'
VERSION = __version__
ISRELEASED = False
########################
CLASSIFIERS = """\
Development Status :: 3 - Alpha
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: MIT License
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Scientific/Engineering :: Bio-Informatics
Topic :: Scientific/Engineering :: Chemistry
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
"""


def find_package_data():
    files = []
    for root, dirnames, filenames in os.walk('pdbfixer'):
        for fn in filenames:
            files.append(relpath(join(root, fn), 'pdbfixer'))
    return files


setup(
    name='pdbfixer',
    author='Peter Eastman',
    description=DOCLINES[0],
    long_description="\n".join(DOCLINES[2:]),
    version=__version__,
    license='MIT',
    url='https://github.com/peastman/pdbfixer',
    platforms=['Linux', 'Mac OS-X', 'Unix', 'Windows'],
    classifiers=CLASSIFIERS.splitlines(),
    packages=find_packages(),
    package_data={'pdbfixer': find_package_data()},
    zip_safe=False,
    entry_points={'console_scripts': ['pdbfixer = pdbfixer.pdbfixer:main']})
    