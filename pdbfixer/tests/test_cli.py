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

import tempfile
import commands

from nose.plugins.skip import Skip, SkipTest

#=============================================================================================
# UNIT TESTS
#=============================================================================================

def run_cli(arguments, expected_output=None):
    [status, output] = commands.getstatusoutput('pdbfixer ' + arguments)

    if status:
        message  = "An error return value (%s) was obtained:\n" % str(status)
        message += "\n"
        message += output
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

def test_url():
    run_cli('--url "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=1LE1"')
