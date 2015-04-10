/**
 * Molecule store.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.store.Molecules', {
  extend: 'Ext.data.Store',
  model: 'Esc.magmaweb.model.Molecule',
  requires: ['Esc.magmaweb.store.proxy.AjaxSingle'],
  proxy: {
    type: 'ajaxsingle',
    listeners: {
      exception: function(proxy, response, operation) {
        if (response.statusText === 'transaction aborted') {
          // a newer request was made so this one was canceled.
          this.fireEvent('abort', proxy, response, operation);
          return;
        }
        Ext.Error.raise({
          msg: 'Failed to load molecules from server',
          response: response,
          operation: operation
        });
      }
    },
    reader: {
      type: 'json',
      root: 'rows',
      idProperty: 'molid'
    }
  },
  sorters: [{
    property: 'refscore',
    direction: 'DESC'
  }, {
    property: 'molid',
    direction: 'ASC'
  }],
  remoteSort: true,
  remoteFilter: true,
  initComponent: function() {
    this.addEvents('abort');
  },
  /**
   * Shortcut for this.getProxy().url
   *
   * @param {String} url URL
   */
  setUrl: function(url) {
    this.getProxy().url = url;
  },
  /**
   * Removes scan filter from molecule store.
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
   * Filter molecules having hits in a specific scan.
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
  setMzFilter: function(mz) {
    this.getProxy().extraParams.mz = mz;
    this.loadPage(1);
  },
  clearMzFilter: function() {
    if ('mz' in this.getProxy().extraParams) {
      delete(this.getProxy().extraParams.mz);
      this.loadPage(1);
    }
  },
  /**
   * Sets the pagesize of the store and reloads store to first page.
   *
   * @param {Number} pageSize Nr of molecules to show on one page.
   */
  setPageSize: function(pageSize) {
    this.pageSize = pageSize;
    this.loadPage(1);
  },
  /**
   * Returns the total number of molecules on server without filtering or paging applied.
   *
   * The `totalUnfiltered` property of the json response.
   *
   * @return {Number} The total number of unfiltered molecules
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
  },
  hasSelectedMolecule: function() {
    return 'molid' in this.getProxy().getReader().rawData;
  }
});