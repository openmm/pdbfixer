import simtk.openmm.app as app
from simtk.openmm.app.internal.pdbstructure import PdbStructure
from pdbfixer import PDBFixer, substitutions
import uiserver
import webbrowser
import os.path
from cStringIO import StringIO

def loadHtmlFile(name):
    htmlPath = os.path.join(os.path.dirname(__file__), 'html')
    file = os.path.join(htmlPath, name)
    return open(file).read()

def startPageCallback(parameters, handler):
    if 'pdbfile' in parameters:
        global fixer
        pdb = PdbStructure(parameters['pdbfile'].value.splitlines())
        fixer = PDBFixer(pdb)
        displayDeleteChainsPage()

def deleteChainsPageCallback(parameters, handler):
    numChains = len(list(fixer.topology.chains()))
    deleteIndices = [i for i in range(numChains) if 'include'+str(i) not in parameters]
    fixer.removeChains(deleteIndices)
    displayAddResiduesPage()

def addResiduesPageCallback(parameters, handler):
    keys = [key for key in sorted(fixer.missingResidues)]
    for i, key in enumerate(keys):
        if 'add'+str(i) not in parameters:
            del fixer.missingResidues[key]
    displayMissingAtomsPage()

def convertResiduesPageCallback(parameters, handler):
    fixer.nonstandardResidues = [residue for i, residue in enumerate(fixer.nonstandardResidues) if 'convert'+str(i) in parameters]
    fixer.replaceNonstandardResidues()
    displayMissingAtomsPage()

def missingAtomsPageCallback(parameters, handler):
    fixer.addMissingAtoms()
    displayAddHydrogensPage()

def addHydrogensPageCallback(parameters, handler):
    if 'add' in parameters:
        pH = float(parameters.getfirst('ph'))
        fixer.addMissingHydrogens(pH)
    displaySaveFilePage()

def saveFilePageCallback(parameters, handler):
    if 'save' in parameters:
        output = StringIO()
        app.PDBFile.writeFile(fixer.topology, fixer.positions, output)
        handler.sendDownload(output.getvalue(), 'output.pdb')
    else:
        displayStartPage()

def displayStartPage():
    uiserver.setCallback(startPageCallback)
    uiserver.setContent(loadHtmlFile("start.html"))

def displayDeleteChainsPage():
    uiserver.setCallback(deleteChainsPageCallback)
    numChains = len(list(fixer.topology.chains()))
    if numChains < 2:
        displayAddResiduesPage()
        return
    table = ""
    proteinResidues = ['ALA', 'ASN', 'CYS', 'GLU', 'HIS', 'LEU', 'MET', 'PRO', 'THR', 'TYR.pdb ARG', 'ASP', 'GLN', 'GLY', 'ILE', 'LYS', 'PHE', 'SER', 'TRP', 'VAL']
    rnaResidues = ['A', 'G', 'C', 'U']
    dnaResidues = ['DA', 'DG', 'DC', 'DT']
    for i, chain in enumerate(fixer.topology.chains()):
        residues = list(r.name for r in chain.residues())
        if any(r in proteinResidues for r in residues):
            content = "Protein"
        elif any(r in rnaResidues for r in residues):
            content = "RNA"
        elif any(r in dnaResidues for r in residues):
            content = "DNA"
        else:
            content = ', '.join(set(residues))
        table += '    <tr><td>%d</td><td>%d</td><td>%s</td><td><input type="checkbox" name="include%d" checked></td></tr>\n' % (chain.index+1, len(residues), content, i)
    uiserver.setContent(loadHtmlFile("removeChains.html") % (numChains, table))

def displayAddResiduesPage():
    uiserver.setCallback(addResiduesPageCallback)
    fixer.findMissingResidues()
    if len(fixer.missingResidues) == 0:
        displayConvertResiduesPage()
        return
    table = ""
    for i, key in enumerate(sorted(fixer.missingResidues)):
        residues = fixer.missingResidues[key]
        table += '    <tr><td>%d</td><td>%d to %d</td><td>%s</td><td><input type="checkbox" name="add%d" checked></td></tr>\n' % (key[0]+1, key[1]+1, key[1]+len(residues), ', '.join(residues), i)
    uiserver.setContent(loadHtmlFile("addResidues.html") % table)

def displayConvertResiduesPage():
    uiserver.setCallback(convertResiduesPageCallback)
    fixer.findNonstandardResidues()
    if len(fixer.nonstandardResidues) == 0:
        displayMissingAtomsPage()
        return
    indexInChain = {}
    for chain in fixer.topology.chains():
        for i, residue in enumerate(chain.residues()):
            indexInChain[residue] = i+1
    table = ""
    for i, residue in enumerate(fixer.nonstandardResidues):
        table += '    <tr><td>%d</td><td>%s %d</td><td>%s</td><td><input type="checkbox" name="convert%d" checked></td></tr>\n' % (residue.chain.index+1, residue.name, indexInChain[residue], substitutions[residue.name], i)
    uiserver.setContent(loadHtmlFile("convertResidues.html") % table)

def displayMissingAtomsPage():
    uiserver.setCallback(missingAtomsPageCallback)
    fixer.findMissingAtoms()
    allResidues = list(set(fixer.missingAtoms.iterkeys()).union(fixer.missingTerminals.iterkeys()))
    allResidues.sort(key=lambda x: x.index)
    if len(allResidues) == 0:
        displayAddHydrogensPage()
        return
    indexInChain = {}
    for chain in fixer.topology.chains():
        for i, residue in enumerate(chain.residues()):
            indexInChain[residue] = i+1
    table = ""
    for residue in allResidues:
        atoms = []
        if residue in fixer.missingAtoms:
            atoms.extend(atom.name for atom in fixer.missingAtoms[residue])
        if residue in fixer.missingTerminals:
            atoms.extend(atom for atom in fixer.missingTerminals[residue])
        table += '    <tr><td>%d</td><td>%s %d</td><td>%s</td></tr>\n' % (residue.chain.index+1, residue.name, indexInChain[residue], ', '.join(atoms))
    uiserver.setContent(loadHtmlFile("addHeavyAtoms.html")% table)

def displayAddHydrogensPage():
    uiserver.setCallback(addHydrogensPageCallback)
    uiserver.setContent(loadHtmlFile("addHydrogens.html"))

def displaySaveFilePage():
    uiserver.setCallback(saveFilePageCallback)
    uiserver.setContent(loadHtmlFile("saveFile.html"))

def launchUI():
    uiserver.beginServing()
    displayStartPage()
    webbrowser.open('http://localhost:'+str(uiserver.server.server_address[1]))