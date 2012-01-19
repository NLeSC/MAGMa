/**
 * Metabolite controller.
 *
 * Handles actions performed in metabolites views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.mmm.controller.Metabolites', {
  extend: 'Ext.app.Controller',
  views: [ 'metabolite.List' ],
  stores: [ 'Metabolites' ],
  models: [ 'Metabolite' ],
  refs: [{
    ref: 'metaboliteList', selector: 'metabolitelist'
  }],
  init: function() {
    console.log('Metabolites controller init');
    var me = this;

    // configure store
    var store = this.getMetabolitesStore();
    store.pageSize = this.application.getPageSize();
    store.setUrl(this.application.getUrls().metabolites);
    store.on('load', this.onLoad, this);

    this.control({
      'metabolitelist': {
        select: this.onSelect,
        deselect: this.onDeselect,
        beforeselect: this.beforeSelect
      },
      'metabolitelist button[action=clear]': {
        click: this.clearFilters
      },
      'metabolitelist component[action=pagesize]': {
        select: this.onPageSizeChange
      },
      'metabolitelist tool[action=download]': {
        click: this.download
      }
    });

    this.application.on('selectscan', this.applyScanFilter, this);
    this.application.on('noselectscan', this.clearScanFilter, this);

    this.application.addEvents(
        /**
         * @event
         * Triggered when metabolite store is loaded.
         * @param {Ext.data.Store} store
         */
        'metaboliteload',
        /**
         * @event
         * Triggered when metabolite is selected.
         * @param {Number} metid Metabolite identifier
         * @param {Esc.mmm.model.Metabolite} metabolite
         */
        'metaboliteselect',
        /**
         * @event
         * Triggered when metabolite is deselected.
         * @param {Number} metid Metabolite identifier
         * @param {Esc.mmm.model.Metabolite} metabolite
         */
        'metabolitedeselect',
        /**
         * @event
         * Triggered when metabolite selection is cleared
         */
        'metabolitenoselect'
    );
  },
  /**
   * Loads metabolite store
   */
  onLaunch: function() {
      // store not loaded in init because metaboliteload event is fired before listeners of views are registerd
      // the nr_scans column has an active filter
      // so do not use list.store.load() , but trigger a filter update to load
      this.getMetaboliteList().filters.createFilters();
      // combo isnt available during init to select pagesize in onlaunch
      Ext.ComponentQuery.query("metabolitelist combo[action=pagesize]")[0].select(this.getMetabolitesStore().pageSize);
  },
  /**
   * Listens for metabolite store load event.
   * Selects metabolite if store only contains 1 metabolite.
   *
   * @param {Ext.data.Store} store
   */
  onLoad: function(store) {
    this.application.fireEvent('metaboliteload', store);
    if (store.getCount() == 1 && !this.getMetaboliteList().getSelectionModel().hasSelection()) {
        console.log('Only one metabolite loaded and its not selected, selecting it');
        this.getMetaboliteList().getSelectionModel().select(0);
    }
  },
  /**
   * Only allow metabolite with a scans to be selected.
   * The extracted ion chromatogram of a metabolite without scans can not be shown because it can not be selected.
   */
  beforeSelect: function(rm, metabolite) {
    return (metabolite.data.nr_scans > 0);
  },
  onSelect: function(rm, metabolite) {
    var metid = metabolite.data.metid;
    this.application.fireEvent('metaboliteselect', metid, metabolite);
  },
  onDeselect: function(rm, metabolite) {
    var metid = metabolite.data.metid;
    this.application.fireEvent('metabolitedeselect', metid, metabolite);
  },
  /**
   * Remove filters and clears selection
   */
  clearFilters: function() {
    console.log('Clear metabolite filters');
    this.getMetaboliteList().clearFilters();
    this.application.fireEvent('metabolitenoselect');
  },
  /**
   * If metabolite is selected then try to reselect it after load
   * @private
   */
  reselectAfterLoad: function() {
      var me = this;
      var sm = me.getMetaboliteList().getSelectionModel();
      var store = me.getMetabolitesStore();
      if (sm.hasSelection()) {
          var selected = sm.getSelection()[0].getId();

          var reselect = function() {
              sm.select(store.getById(selected));
              store.removeListener('load', reselect, me);
          }

          store.on('load', reselect , me);
      }
  },
  /**
   * Apply scan filter to metabolite store.
   * And shows fragment score column.
   * Tries to keep selection.
   *
   * @param {Number} scanid Scan identifier to filter on.
   */
  applyScanFilter: function(scanid) {
      this.reselectAfterLoad();
      this.getMetabolitesStore().setScanFilter(scanid);
      this.getMetaboliteList().getFragmentScoreColumn().show();
      // TODO in spectra add markers for metabolites present in scan
  },
  /**
   * Removes scan filter from metabolite store.
   * And hides fragment score column.
   * And deactivates filters on fragment score column if any
   * And resets sort if store is sorted on fragment score.
   */
  clearScanFilter: function() {
      var store = this.getMetabolitesStore();
      if ('filters' in store) {
          store.filters.removeAtKey('score');
      }
      if ('sorters' in store) {
          store.sorters.removeAtKey('score');
      }
      store.removeScanFilter();
      this.getMetaboliteList().getFragmentScoreColumn().hide();
  },
  onPageSizeChange: function(combo) {
      this.getMetabolitesStore().setPageSize(combo.getValue());
  },
  /**
   * Open a new window with metabolites as comma seperated file.
   * Uses application.urls.metabolitescsv as url.
   * Uses store/proxy/gridfilter state to construct queryString so what you see in grid is what in csv file.
   */
  download: function() {
    // download needs to make an url with has the same query parameters as the store.load()
    // for load() the store builds an operation object, we need to build this aswell
    var store = this.getMetabolitesStore();
    var proxy = store.getProxy();
    var config = {
        action: 'read',
        groupers: store.groupers.items,
        limit: store.pageSize,
        page: store.currentPage,
        start: (store.currentPage-1)*store.pageSize,
        sorters: store.getSorters()
    };
    var operation = Ext.create('Ext.data.Operation', config);
    var request = proxy.buildRequest(operation);
    var params = request.params;
    // Ext.ux.grid.FiltersFeature adds filters to request.params in store.beforeLoad event handler
    // so we do the same to get the filter query string
    var filter = this.getMetaboliteList().getView().getFeature('mfilter');
    Ext.apply(params, filter.buildQuery(filter.getFilterData()));

    var url = Ext.urlAppend(
        this.application.getUrls().metabolitescsv,
        Ext.Object.toQueryString(params)
    );
    window.open(url, 'metabolites.csv');
  }
});

