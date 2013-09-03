import simtk.openmm.app as app
from pdbfixer import PDBFixer, substitutions
import uiserver
import webbrowser
from cStringIO import StringIO

def startPageCallback(parameters, handler):
    if 'pdbfile' in parameters:
        global fixer
        pdb = app.PDBFile(parameters['pdbfile'].value.splitlines())
        fixer = PDBFixer(pdb)
        displayConvertResiduesPage()

def convertResiduesPageCallback(parameters, handler):
    global nonstandard
    fixer.replaceNonstandardResidues(nonstandard)
    nonstandard = []
    displayMissingAtomsPage()

def missingAtomsPageCallback(parameters, handler):
    fixer.addMissingAtoms(missingAtoms, missingTerminals)
    displayDownloadPage()

def downloadPageCallback(parameters, handler):
    if 'download' in parameters:
        output = StringIO()
        app.PDBFile.writeFile(fixer.processedTopology, fixer.processedPositions, output)
        handler.sendDownload(output.getvalue(), 'output.pdb')
    else:
        displayStartPage()

def displayStartPage():
    uiserver.setCallback(startPageCallback)
    uiserver.setContent("""
<html>
<head><title>PDBFixer</title></head>
<body>
<h1>Welcome To PDBFixer!</h1>
Select a PDB file to load.  It will be analyzed for problems.
<p>
<form method="post" action="/" enctype="multipart/form-data">
PDB File: <input type="file" name="pdbfile"/>
<p>
<input type="submit" value="Analyze File"/>
</form>
</body>
</html>
""")

def displayConvertResiduesPage():
    uiserver.setCallback(convertResiduesPageCallback)
    global nonstandard
    nonstandard = fixer.findNonstandardResidues()
    if len(nonstandard) == 0:
        displayMissingAtomsPage()
        return
    table = ""
    for i, residue in enumerate(nonstandard):
        table += '    <tr><td>%d</td><td>%s %d</td><td>%s</td><td><input type="checkbox" name="convert%d" checked></td></tr>\n' % (residue.chain.index+1, residue.name, residue.index+1, substitutions[residue.name], i)
    uiserver.setContent("""
<html>
<head><title>PDB Fixer</title></head>
<body>
This PDB file contains non-standard residues.  Do you want to convert them to the corresponding standard residues?
<p>
<form method="get" action="/">
<table border="1">
    <tr><th>Chain</th><th>Residue</th><th>Convert To</th><th>Convert?</th></tr>
%s
</table>
<p>
<input type="submit" value="Continue"/>
</form>
</body>
<html>
""" % table)

def displayMissingAtomsPage():
    uiserver.setCallback(missingAtomsPageCallback)
    global missingAtoms
    global missingTerminals
    missingAtoms, missingTerminals = fixer.findMissingAtoms()
    allResidues = list(set(missingAtoms.iterkeys()).union(missingTerminals.iterkeys()))
    allResidues.sort(key=lambda x: x.index)
    if len(allResidues) == 0:
        displayDownloadPage()
        return
    table = ""
    for residue in allResidues:
        atoms = []
        if residue in missingAtoms:
            atoms.extend(atom.name for atom in missingAtoms[residue])
        if residue in missingTerminals:
            atoms.extend(atom for atom in missingTerminals[residue])
        table += '    <tr><td>%d</td><td>%s %d</td><td>%s</td></tr>\n' % (residue.chain.index+1, residue.name, residue.index+1, ', '.join(atoms))
    uiserver.setContent("""
<html>
<head><title>PDB Fixer</title></head>
<body>
The following residues are missing heavy atoms, which will be added.
<p>
<form method="get" action="/">
<table border="1">
    <tr><th>Chain</th><th>Residue</th><th>Missing Atoms</th></tr>
%s
</table>
<p>
<input type="submit" value="Continue"/>
</form>
</body>
<html>
""" % table)

def displayDownloadPage():
    uiserver.setCallback(downloadPageCallback)
    uiserver.setContent("""
<html>
<head><title>PDB Fixer</title></head>
<body>
The fixed PDB file is ready to download.
<p>
<form method="get" action="/">
<input type="submit" name="download" value="Download"/>
<input type="submit" name="newfile" value="Process Another File"/>
</form>
</body>
<html>
""")

def launchUI():
    uiserver.beginServing()
    displayStartPage()
    webbrowser.open('http://'+uiserver.serverAddress)