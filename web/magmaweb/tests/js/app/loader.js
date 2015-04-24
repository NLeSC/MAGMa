// disable dynamic dependency loading, karma loads all the required files in the correct order
Ext.Loader.setConfig({
  enabled: false
});
// running in test environment somehow makes Ext.resetElement undefined.
if (!('resetElement' in Ext)) {
  Ext.resetElement = Ext.getBody();
}

// setup which normally resides in index.html
Ext.getBody().createChild({
  tag: 'div',
  id: 'resultsinfo',
  html: 'placeholder for resultsinfo'
});
Ext.getBody().createChild({
  tag: 'div',
  id: 'logos',
  html: 'placeholder for logos'
});

var sketcher = {getMolecule: function() { return {bonds:[]}; } };

var appRootBase = 'base/magmaweb/tests/js/app';