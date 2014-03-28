var glmol = null;
var defaultRepresentation = {
    color: 'chainbow',
    mainChain: 'thickRibbon',
    sideChains: 'hide',
    nonbondedAtoms: 'stars',
    heteroAtoms: 'sphere',
    backgroundColor: 0xFFFFFF,
}

renderPDB = function(pdbdata) {
    glmol.loadMoleculeStr(false, pdbdata);
    glmol.rebuildScene(defaultRepresentation);
    glmol.show();
};

$(function() {
    container = $('#proteinViewContainer');
    if ($('#proteinViewContainer').length > 0) {
        $('#mainform').css({float: 'left'});
        container.height(500);
        container.width(500);
        container.css({float: 'left'});
        container.css({'border': '1px black'})
        $('#proteinViewContainer').width(500);
        
        $('<a>').attr('href', '#')
        .html('View Molecule')
        .click(function(){
            if (glmol == null) {
                glmol = new GLmol('proteinView', true);
            }
            $.get('/current-pdb', function(data) {
                renderPDB(data);
            });
        }).appendTo($('body'));
    }
    
});