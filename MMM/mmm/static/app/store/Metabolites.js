/**
 * Metabolite store.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.mmm.store.Metabolites', {
  extend: 'Ext.data.Store',
  model: 'Esc.mmm.model.Metabolite',
  proxy: {
    type: 'ajax',
    reader: {
      type: 'json',
      root: 'rows',
      idProperty: 'metid'
    }
  },
  sorters: [{
    property: 'probability',
    direction: 'DESC'
  },{
    property: 'metid',
    direction: 'ASC'
  }],
  remoteSort: true,
  remoteFilter: true,
  isLoaded: false,
  /**
   * Shortcut for this.getProxy().url
   *
   * @param {String} url URL
   */
  setUrl: function(url) {
    this.getProxy().url = url;
  },
  listeners: {
    load: function(store) {
      this.isLoaded = true;
    }
  },
  /**
   * Removes scan filter from metabolite store.
   * And reloads store to first page.
   */
  removeScanFilter: function() {
    // see if already filtered on scanid then remove old filter
    if ('scanid' in this.getProxy().extraParams) {
      delete(this.getProxy().extraParams.scanid);
      this.loadPage(1);
    }
  },
  /**
   * Filter metabolites having hits in a specific scan.
   * Sets filter and reloads store to first page.
   *
   * @param {Number} scanid Scan identifier.
   */
  setScanFilter: function(scanid) {
    this.getProxy().extraParams.scanid = scanid;
    this.loadPage(1);
  }
});
