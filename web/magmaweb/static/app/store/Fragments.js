/**
 * Store for fragments.
 *
 * Fragments are loaded when a scan and molecule are selected.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.store.Fragments', {
  extend: 'Ext.data.TreeStore',
  model: 'Esc.magmaweb.model.Fragment',
  autoLoad: false,
  clearOnLoad: false,
  root: {
    expanded: false,
    children: []
  }, // prevent tree from autoloading
  /**
   * TreeStore and Store have different function to fetch record by id, add getById to TreeStore
   *
   * @param {Number} id Identifier of fragment
   * @return {Ext.data.NodeInterface}
   */
  getById: function(id) {
    return this.getNodeById(id);
  },
  /**
   * Find a fragment by m/z and MS level.
   *
   * @param {Number} mz M/z of node fragment to find
   * @param {Number} mslevel MS level on which m/z must be found
   * @return {Ext.data.NodeInterface}
   */
  getNodeByMzMslevel: function(mz, mslevel) {
    return this.getRootNode().findChildBy(function(n) {
      return (n.data.mslevel == mslevel && n.data.mz == mz);
    }, false, true);
  }
});