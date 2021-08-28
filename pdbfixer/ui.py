from __future__ import absolute_import

import webbrowser
import os.path
import time
import sys
from math import sqrt

import openmm.app as app
import openmm.unit as unit
from openmm.vec3 import Vec3

from .pdbfixer import PDBFixer, proteinResidues, dnaResidues, rnaResidues, _guessFileFormat
from . import uiserver

if sys.version_info >= (3,0):
    from io import StringIO
else:
    from cStringIO import StringIO

def loadHtmlFile(name):
    htmlPath = os.path.join(os.path.dirname(__file__), 'html')
    file = os.path.join(htmlPath, name)
    return open(file).read()

cachedImages = {}

def loadImageFile(name):
    global cachedImages
    if name not in cachedImages:
        imagePath = os.path.join(os.path.dirname(__file__), 'images')
        file = os.path.join(imagePath, name)
        cachedImages[name] = open(file, 'rb').read()
    return cachedImages[name]

def controlsCallback(parameters, handler):
    if 'newfile' in parameters:
        displayStartPage()
    if 'quit' in parameters:
        handler.sendResponse(loadHtmlFile("quit.html"))
        uiserver.server.shutdown()
        global uiIsRunning
        uiIsRunning = False

def imageCallback(parameters, handler):
    name = parameters['name'][0]
    image = loadImageFile(name)
    type = None
    if name.endswith('.png'):
        type = 'image/png'
    elif name.endswith('.jpeg') or name.endswith('.jpg'):
        type = 'image/jpeg'
    handler.sendResponse(image, type=type)

def startPageCallback(parameters, handler):
    global fixer
    if 'type' in parameters:
        if parameters.getfirst('type') == 'local':
            filename = parameters['pdbfile'].filename
            file = StringIO(parameters['pdbfile'].value.decode())
            if _guessFileFormat(file, filename) == 'pdbx':
                fixer = PDBFixer(pdbxfile=file)
            else:
                fixer = PDBFixer(pdbfile=file)
            fixer.source = filename
        else:
            id = parameters.getfirst('pdbid')
            try:
                fixer = PDBFixer(pdbid=id)
            except Exception as e:
                import traceback
                print(traceback.format_exc())
                handler.sendResponse(
                    header + "<p>Unable to download the PDB file. " +
                    "This may indicate an invalid PDB identifier, " +
                    "or an error in network connectivity.</p>" +
                    "<p>{}</p>".format(e) +
                    loadHtmlFile("error.html"))
        displayDeleteChainsPage()

def deleteChainsPageCallback(parameters, handler):
    global heterogens
    heterogens = parameters.getfirst('heterogens')
    numChains = len(list(fixer.topology.chains()))
    deleteIndices = [i for i in range(numChains) if 'include'+str(i) not in parameters]
    fixer.removeChains(deleteIndices)
    displayAddResiduesPage()

def addResiduesPageCallback(parameters, handler):
    keys = [key for key in sorted(fixer.missingResidues)]
    for i, key in enumerate(keys):
        if 'add'+str(i) not in parameters:
            del fixer.missingResidues[key]
    displayConvertResiduesPage()

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
    if 'addhydrogens' in parameters:
        pH = float(parameters.getfirst('ph'))
        fixer.addMissingHydrogens(pH)
    if 'addwater' in parameters:
        padding, boxSize, boxVectors = None, None, None
        if parameters.getfirst('boxType') == 'geometry':
            geompadding = float(parameters.getfirst('geomPadding')) * unit.nanometer
            geometry = parameters.getfirst('geometryDropdown')
            base_size = float(parameters.getfirst('maxMolecularAxis')) * unit.nanometer
            if geometry == 'cube':
                padding = geompadding
            elif geometry == 'truncatedOctahedron':
                vectors = Vec3(1,0,0), Vec3(1/3,2*sqrt(2)/3,0), Vec3(-1/3,sqrt(2)/3,sqrt(6)/3)
                boxVectors = [(base_size+geompadding)*v for v in vectors]
            elif geometry == 'rhombicDodecahedron':
                vectors = Vec3(1,0,0), Vec3(0,1,0), Vec3(0.5,0.5,sqrt(2)/2)
                boxVectors = [(base_size+geompadding)*v for v in vectors]
        else:
            boxSize = (float(parameters.getfirst('boxx')), float(parameters.getfirst('boxy')), float(parameters.getfirst('boxz')))*unit.nanometer
        ionicStrength = float(parameters.getfirst('ionicstrength'))*unit.molar
        positiveIon = parameters.getfirst('positiveion')+'+'
        negativeIon = parameters.getfirst('negativeion')+'-'
        fixer.addSolvent(boxSize, padding, boxVectors, positiveIon, negativeIon, ionicStrength)
    if 'addmembrane' in parameters:
        lipidType = parameters.getfirst('lipidType')
        padding = float(parameters.getfirst('membranePadding'))*unit.nanometer
        ionicStrength = float(parameters.getfirst('ionicstrength'))*unit.molar
        positiveIon = parameters.getfirst('positiveion')+'+'
        negativeIon = parameters.getfirst('negativeion')+'-'
        fixer.addMembrane(lipidType, 0*unit.nanometer, padding, positiveIon, negativeIon, ionicStrength)
    displaySaveFilePage()

def saveFilePageCallback(parameters, handler):
    if 'pdb' in parameters:
        output = StringIO()
        if fixer.source is not None:
            output.write("REMARK   1 PDBFIXER FROM: %s\n" % fixer.source)
        try:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, output, True)
        except AssertionError:
            print
        handler.sendDownload(output.getvalue(), 'output.pdb')
    elif 'pdbx' in parameters:
        output = StringIO()
        if fixer.source is not None:
            output.write("# Created with PDBFixer from: %s\n" % fixer.source)
        try:
            app.PDBxFile.writeFile(fixer.topology, fixer.positions, output, True)
        except AssertionError:
            print
        handler.sendDownload(output.getvalue(), 'output.cif')
    else:
        displayStartPage()

def displayStartPage():
    uiserver.setCallback(startPageCallback)
    uiserver.setContent(header+loadHtmlFile("start.html"))

def displayDeleteChainsPage():
    uiserver.setCallback(deleteChainsPageCallback)
    numChains = len(list(fixer.topology.chains()))
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
        table += '    <tr><td>%s</td><td>%d</td><td>%s</td><td><input type="checkbox" name="include%d" checked></td></tr>\n' % (chain.id, len(residues), content, i)
    uiserver.setContent(header+loadHtmlFile("removeChains.html") % (numChains, table))

def displayAddResiduesPage():
    uiserver.setCallback(addResiduesPageCallback)
    fixer.findMissingResidues()
    if len(fixer.missingResidues) == 0:
        displayConvertResiduesPage()
        return
    table = ""
    chains = list(fixer.topology.chains())
    for i, key in enumerate(sorted(fixer.missingResidues)):
        residues = fixer.missingResidues[key]
        chain = chains[key[0]]
        chainResidues = list(chain.residues())
        if key[1] < len(chainResidues):
            offset = int(chainResidues[key[1]].id)-len(residues)-1
        else:
            offset = int(chainResidues[-1].id)
        table += '    <tr><td>%s</td><td>%d to %d</td><td>%s</td><td><input type="checkbox" name="add%d" checked></td></tr>\n' % (chain.id, offset+1, offset+len(residues), ', '.join(residues), i)
    uiserver.setContent(header+loadHtmlFile("addResidues.html") % table)

def displayConvertResiduesPage():
    uiserver.setCallback(convertResiduesPageCallback)
    fixer.findNonstandardResidues()
    if len(fixer.nonstandardResidues) == 0:
        displayMissingAtomsPage()
        return
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
        table += '    <tr><td>%s</td><td>%s %s</td><td><select name="residue%d">%s</select></td><td><input type="checkbox" name="convert%d" checked></td></tr>\n' % (residue.chain.id, residue.name, residue.id, i, options, i)
    uiserver.setContent(header+loadHtmlFile("convertResidues.html") % table)

def displayMissingAtomsPage():
    uiserver.setCallback(missingAtomsPageCallback)
    if heterogens == 'none':
        fixer.removeHeterogens(False)
    elif heterogens == 'water':
        fixer.removeHeterogens(True)
    fixer.findMissingAtoms()
    allResidues = list(set(fixer.missingAtoms.keys()).union(fixer.missingTerminals.keys()))
    allResidues.sort(key=lambda x: x.index)
    if len(allResidues) == 0:
        fixer.addMissingAtoms()
        displayAddHydrogensPage()
        return
    table = ""
    for residue in allResidues:
        atoms = []
        if residue in fixer.missingAtoms:
            atoms.extend(atom.name for atom in fixer.missingAtoms[residue])
        if residue in fixer.missingTerminals:
            atoms.extend(atom for atom in fixer.missingTerminals[residue])
        table += '    <tr><td>%s</td><td>%s %s</td><td>%s</td></tr>\n' % (residue.chain.id, residue.name, residue.id, ', '.join(atoms))
    uiserver.setContent(header+loadHtmlFile("addHeavyAtoms.html") % table)

def displayAddHydrogensPage():
    uiserver.setCallback(addHydrogensPageCallback)
    dimensions = ""
    if fixer.topology.getUnitCellDimensions() is not None:
        dimensions = "<tr><td>Crystallographic unit cell:</td><td>%.3f</td><td>%.3f</td><td>%.3f</td></tr>" % fixer.topology.getUnitCellDimensions().value_in_unit(unit.nanometer)
    sizeRange = tuple(max((pos[i] for pos in fixer.positions))-min((pos[i] for pos in fixer.positions)) for i in range(3))
    dimensions += "<tr id='boxContainingAllAtoms'><td>Box containing all atoms:</td><td>%.3f</td><td>%.3f</td><td>%.3f</td></tr>" % tuple(x.value_in_unit(unit.nanometer) for x in sizeRange)
    uiserver.setContent(header+loadHtmlFile("addHydrogens.html") % dimensions)

def displaySaveFilePage():
    uiserver.setCallback(saveFilePageCallback)
    uiserver.setContent(header+loadHtmlFile("saveFile.html"))

def launchUI():
    global header
    header = loadHtmlFile("header.html")
    uiserver.beginServing()
    uiserver.setCallback(controlsCallback, "/controls")
    uiserver.setCallback(imageCallback, "/image")
    displayStartPage()
    url = 'http://localhost:'+str(uiserver.server.server_address[1])
    print("PDBFixer running: %s " % url)
    webbrowser.open(url)

    # the uiserver is running in a background daemon thread that dies whenever
    # the main thread exits. So, to keep the whole process alive, we just sleep
    # here in the main thread. When Control-C is called, the main thread shuts
    # down and then the uiserver exits. Without this daemon/sleep combo, the
    # process cannot be killed with Control-C. Reference stack overflow link:
    # http://stackoverflow.com/a/11816038/1079728

    global uiIsRunning
    uiIsRunning = True
    while uiIsRunning:
        time.sleep(0.5)
