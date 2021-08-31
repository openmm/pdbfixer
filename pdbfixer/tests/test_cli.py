#!/usr/bin/python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Test command-line interface.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import subprocess
from subprocess import CalledProcessError

import pytest

#=============================================================================================
# UNIT TESTS
#=============================================================================================

def run_cli(arguments, expected_output=None):
    try:
        output = subprocess.check_output('pdbfixer ' + arguments, shell=True)
    except CalledProcessError as e:
        message  = "An error return value (%s) was obtained:\n" % str(e.returncode)
        message += "\n"
        message += str(e.output)
        message += "\n"
        raise Exception(message)

    if expected_output:
        if output != expected_output:
            message  = "Output differs from expected output.\n"
            message += "\n"
            message += "Expected output:\n"
            message += expected_output
            message += "\n"
            message += "Actual output:\n"
            message += output
            message += "\n"
            raise Exception(message)

def test_help():
    run_cli('--help')

def test_pdbid():
    run_cli('--pdbid 1LE1')

@pytest.mark.skipif(os.getenv("GITHUB_ACTION") is not None, reason="Cannot download during CI")
def test_url():
    run_cli('--url "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=1LE1"')
