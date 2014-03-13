import simtk.openmm.app as app
from simtk.openmm.app.internal.pdbstructure import PdbStructure
import simtk.unit as unit
from pdbfixer import PDBFixer, substitutions, proteinResidues, dnaResidues, rnaResidues
import uiserver
import webbrowser
import os.path
import gzip
from io import BytesIO
try:
    from urllib.request import urlopen
    from io import StringIO
except:
    from urllib2 import urlopen
    from cStringIO import StringIO

def loadHtmlFile(name):
    htmlPath = os.path.join(os.path.dirname(__file__), 'html')
    file = os.path.join(htmlPath, name)
    return open(file).read()

def controlsCallback(parameters, handler):
    if 'newfile' in parameters:
        displayStartPage()
    if 'quit' in parameters:
        handler.sendResponse(loadHtmlFile("quit.html"))
        uiserver.server.shutdown()

def startPageCallback(parameters, handler):
    global fixer
    if 'type' in parameters:
        if parameters.getfirst('type') == 'local':
            pdb = PdbStructure(parameters['pdbfile'].value.decode().splitlines())
            fixer = PDBFixer(pdb)
            displayDeleteChainsPage()
        else:
            id = parameters.getfirst('pdbid')
            url = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb"+id.lower()+".ent.gz"
            try:
                response = urlopen(url)
                content = gzip.GzipFile(fileobj=BytesIO(response.read())).read()
                pdb = PdbStructure(content.decode().splitlines())
                fixer = PDBFixer(pdb)
                displayDeleteChainsPage()
            except:
                handler.sendResponse(header+"Unable to download the PDB file. This may indicate an invalid PDB identifier, or an error in network connectivity."+loadHtmlFile("error.html"))

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
    for i in range(len(fixer.nonstandardResidues)):
        if 'convert'+str(i) in parameters:
            fixer.nonstandardResidues[i] = (fixer.nonstandardResidues[i][0], parameters.getfirst('residue'+str(i)))
    fixer.replaceNonstandardResidues()
    displayMissingAtomsPage()

def missingAtomsPageCallback(parameters, handler):
    fixer.addMissingAtoms()
    displayAddHydrogensPage()

def addHydrogensPageCallback(parameters, handler):
    heterogens = parameters.getfirst('heterogens')
    if heterogens == 'none':
        fixer.removeHeterogens(False)
    elif heterogens == 'water':
        fixer.removeHeterogens(True)
    if 'addhydrogens' in parameters:
        pH = float(parameters.getfirst('ph'))
        fixer.addMissingHydrogens(pH)
    if 'addwater' in parameters:
        boxSize = (float(parameters.getfirst('boxx')), float(parameters.getfirst('boxy')), float(parameters.getfirst('boxz')))*unit.nanometer
        ionicStrength = float(parameters.getfirst('ionicstrength'))*unit.molar
        positiveIon = parameters.getfirst('positiveion')+'+'
        negativeIon = parameters.getfirst('negativeion')+'-'
        fixer.addSolvent(boxSize, positiveIon, negativeIon, ionicStrength)
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
    uiserver.setContent(header+loadHtmlFile("start.html"))

def displayDeleteChainsPage():
    uiserver.setCallback(deleteChainsPageCallback)
    numChains = len(list(fixer.topology.chains()))
    if numChains < 2:
        displayAddResiduesPage()
        return
    table = ""
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
    uiserver.setContent(header+loadHtmlFile("removeChains.html") % (numChains, table))

def displayAddResiduesPage():
    uiserver.setCallback(addResiduesPageCallback)
    fixer.findMissingResidues()
    if len(fixer.missingResidues) == 0:
        displayConvertResiduesPage()
        return
    table = ""
    for i, key in enumerate(sorted(fixer.missingResidues)):
        residues = fixer.missingResidues[key]
        chain = fixer.structureChains[key[0]]
        if key[1] < len(chain.residues):
            offset = chain.residues[key[1]].number-len(residues)-1
        else:
            offset = chain.residues[-1].number
        table += '    <tr><td>%d</td><td>%d to %d</td><td>%s</td><td><input type="checkbox" name="add%d" checked></td></tr>\n' % (key[0]+1, offset+1, offset+len(residues), ', '.join(residues), i)
    uiserver.setContent(header+loadHtmlFile("addResidues.html") % table)

def displayConvertResiduesPage():
    uiserver.setCallback(convertResiduesPageCallback)
    fixer.findNonstandardResidues()
    if len(fixer.nonstandardResidues) == 0:
        displayMissingAtomsPage()
        return
    indexInChain = {}
    for structChain, topChain in zip(fixer.structureChains, fixer.topology.chains()):
        for structResidue, topResidue in zip(structChain.iter_residues(), topChain.residues()):
            indexInChain[topResidue] = structResidue.number
    table = ''
    nucleotides = ['DA', 'DC', 'DG', 'DT', 'A', 'C', 'G', 'T']
    for i in range(len(fixer.nonstandardResidues)):
        residue, replaceWith = fixer.nonstandardResidues[i]
        if replaceWith in proteinResidues:
            replacements = proteinResidues
        else:
            replacements = nucleotides
        options = ''
        for res in replacements:
            selected = ''
            if res == replaceWith:
                selected = ' selected'
            options += '<option value="%s"%s>%s</option>' % (res, selected, res)
        table += '    <tr><td>%d</td><td>%s %d</td><td><select name="residue%d">%s</select></td><td><input type="checkbox" name="convert%d" checked></td></tr>\n' % (residue.chain.index+1, residue.name, indexInChain[residue], i, options, i)
    uiserver.setContent(header+loadHtmlFile("convertResidues.html") % table)

def displayMissingAtomsPage():
    uiserver.setCallback(missingAtomsPageCallback)
    fixer.findMissingAtoms()
    allResidues = list(set(fixer.missingAtoms.keys()).union(fixer.missingTerminals.keys()))
    allResidues.sort(key=lambda x: x.index)
    if len(allResidues) == 0:
        fixer.addMissingAtoms()
        displayAddHydrogensPage()
        return
    indexInChain = {}
    for structChain, topChain in zip(fixer.structureChains, fixer.topology.chains()):
        for structResidue, topResidue in zip(structChain.iter_residues(), topChain.residues()):
            indexInChain[topResidue] = structResidue.number
    table = ""
    for residue in allResidues:
        atoms = []
        if residue in fixer.missingAtoms:
            atoms.extend(atom.name for atom in fixer.missingAtoms[residue])
        if residue in fixer.missingTerminals:
            atoms.extend(atom for atom in fixer.missingTerminals[residue])
        table += '    <tr><td>%d</td><td>%s %d</td><td>%s</td></tr>\n' % (residue.chain.index+1, residue.name, indexInChain[residue], ', '.join(atoms))
    uiserver.setContent(header+loadHtmlFile("addHeavyAtoms.html") % table)

def displayAddHydrogensPage():
    uiserver.setCallback(addHydrogensPageCallback)
    dimensions = ""
    if fixer.topology.getUnitCellDimensions() is not None:
        dimensions = "<tr><td>Crystallographic unit cell:</td><td>%.3f</td><td>%.3f</td><td>%.3f</td></tr>" % fixer.topology.getUnitCellDimensions().value_in_unit(unit.nanometer)
    sizeRange = tuple(max((pos[i] for pos in fixer.positions))-min((pos[i] for pos in fixer.positions)) for i in range(3))
    dimensions += "<tr><td>Box containing all atoms:</td><td>%.3f</td><td>%.3f</td><td>%.3f</td></tr>" % tuple(x.value_in_unit(unit.nanometer) for x in sizeRange)
    uiserver.setContent(header+loadHtmlFile("addHydrogens.html") % dimensions)

def displaySaveFilePage():
    uiserver.setCallback(saveFilePageCallback)
    uiserver.setContent(header+loadHtmlFile("saveFile.html"))

def launchUI():
    global header
    header = loadHtmlFile("header.html")
    uiserver.beginServing()
    uiserver.setCallback(controlsCallback, "/controls")
    displayStartPage()
    webbrowser.open('http://localhost:'+str(uiserver.server.server_address[1]))
