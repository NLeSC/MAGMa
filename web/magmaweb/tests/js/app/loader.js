// disable dynamic dependency loading, karma loads all the required files in the correct order
Ext.Loader.setConfig({
  enabled: false
});
// running in test environment somehow makes Ext.resetElement undefined.
if (!('resetElement' in Ext)) {
  Ext.resetElement = Ext.getBody();
}