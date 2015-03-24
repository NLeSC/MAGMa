/**
 * Metabolite store.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.store.Metabolites', {
  extend: 'Ext.data.Store',
  model: 'Esc.magmaweb.model.Metabolite',
  proxy: {
    type: 'ajax',
    listeners: {
      exception: function(proxy, response, operation) {
          Ext.Error.raise({
              msg: 'Failed to load metabolites from server',
              response: response,
              operation: operation
          });
      }
    },
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
  /**
   * Shortcut for this.getProxy().url
   *
   * @param {String} url URL
   */
  setUrl: function(url) {
    this.getProxy().url = url;
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
  },
  isFilteredOnScan: function() {
    return 'scanid' in this.getProxy().extraParams;
  },
  /**
   * Sets the pagesize of the store and reloads store to first page.
   *
   * @param {Number} pageSize Nr of metabolites to show on one page.
   */
  setPageSize: function(pageSize) {
    this.pageSize = pageSize;
    this.loadPage(1);
  },
  /**
   * Returns the total number of metabolites on server without filtering or paging applied.
   *
   * The `totalUnfiltered` property of the json response.
   *
   * @return {Number} The total number of unfiltered metabolites
   */
  getTotalUnfilteredCount: function() {
      var reader = this.getProxy().getReader();
      if ('rawData' in reader) {
          return reader.rawData.totalUnfiltered || 0;
      } else {
          return 0;
      }
  },
  selectMolecule: function(molecule) {
      this.getProxy().extraParams.molid = molecule.getId();
  },
  clearMoleculeSelection: function() {
      delete(this.getProxy().extraParams.molid);
  }
});
