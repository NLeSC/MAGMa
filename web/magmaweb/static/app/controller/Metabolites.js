/**
 * Metabolite controller.
 *
 * Handles actions performed in metabolites views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.controller.Metabolites', {
  extend: 'Ext.app.Controller',
  views: [ 'metabolite.Panel' ],
  stores: [ 'Metabolites' ],
  models: [ 'Metabolite' ],
  uses: [
    'Esc.magmaweb.view.metabolite.MetabolizeForm',
    'Esc.magmaweb.view.metabolite.MetabolizeOneForm'
  ],
  refs: [{
    ref: 'metaboliteList', selector: 'metabolitelist'
  }, {
    ref: 'metaboliteAddForm', selector: 'metaboliteaddform'
  }, {
    ref: 'metabolitePanel', selector: 'metabolitepanel'
  }],
  init: function() {
    Ext.log({}, 'Metabolites controller init');
    var me = this;

    // configure store
    var store = this.getMetabolitesStore();
    store.pageSize = this.application.getPageSize();
    store.setUrl(this.application.metabolitesUrl('json'));
    store.on('load', this.onLoad, this);
    store.on('beforeload', this.onBeforeLoad, this);

    this.control({
      'metabolitelist': {
        select: this.onSelect,
        deselect: this.onDeselect,
        beforeselect: this.beforeSelect,
        metabolize: this.showMetabolizeStructureForm
      },
      'metabolitelist component[action=pagesizeCombo]': {
        select: this.onPageSizeChange
      },
      'metabolitepanel component[action=download]': {
        click: this.showDownloadMenu
      },
      'metabolitepanel component[action=actions]': {
        click: this.showActionsMenu
      },
      'metaboliteaddform component[action=addstructures]': {
        click: this.addStructuresHandler
      },
      'metaboliteaddform component[action=addstructurescancel]': {
        click: this.showGrid
      },
      'metabolizeform component[action=metabolize]': {
        click: this.metabolizeHandler
      },
      'metabolizeoneform component[action=metabolize_one]': {
        click: this.metabolizeOneHandler
      }
    });

    this.application.on('selectscan', this.applyScanFilter, this);
    this.application.on('noselectscan', this.clearScanFilter, this);

    /**
     * @property {Boolean} hasMSData
     * Whether there is ms data to use for annotate
     * Used in action forms to disable/enable annotate options
     */
    this.hasMSData = false;
    this.application.on('chromatogramload', this.onChromatrogramLoad, this);
    this.application.on('rpcsubmitsuccess', function() {
      Ext.getCmp('addstructuresaction').disable();
      Ext.getCmp('metabolizeaction').disable();
      me.getMetaboliteList().getCommandsColumn().disableAction();
    });

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
         * @param {Esc.magmaweb.model.Metabolite} metabolite
         */
        'metaboliteselect',
        /**
         * @event
         * Triggered when metabolite is deselected.
         * @param {Number} metid Metabolite identifier
         * @param {Esc.magmaweb.model.Metabolite} metabolite
         */
        'metabolitedeselect',
        /**
         * @event
         * Triggered when metabolite selection is cleared
         */
        'metabolitenoselect'
    );

    this.actionsMenu = Ext.create('Ext.menu.Menu', {
        items: [{
            iconCls: 'icon-add',
            id: 'addstructuresaction',
            text: 'Add structures',
            handler: this.showAddStructuresForm.bind(this)
        }, {
            text: 'Metabolize',
            id: 'metabolizeaction',
            tooltip: 'Metabolize all structures',
            disabled: true,
            handler: this.showMetabolizeForm.bind(this)
        }, {
            text: 'Clear filters',
            handler: this.clearFilters.bind(this)
        }]
    });

    this.downloadMenu = Ext.create('Ext.menu.Menu', {
        items: [{
            text: 'CSV',
            tooltip: 'Save metabolites as comma seperated file',
            handler: this.download_csv.bind(this)
        }, {
            text: 'SDF',
            tooltip: 'Save metabolites as sdf',
            handler: this.download_sdf.bind(this)
        }]
    });

    this.application.on('assignmentchanged', function(isAssigned, params) {
        me.getMetabolitesStore().load();
    });
  },
  /**
   * Loads metabolite store
   */
  onLaunch: function() {
      // store not loaded in init because metaboliteload event is fired before listeners of views are registerd
      // the nhits column has an active filter
      // so do not use list.store.load() , but trigger a filter update to load
      this.getMetaboliteList().filters.createFilters();
      this.getMetabolitesStore().load();
      this.getMetaboliteList().setPageSize(this.getMetabolitesStore().pageSize);
      this.applyRole();
  },
  /**
   * Listens for metabolite store load event.
   * Selects metabolite if store only contains 1 metabolite.
   *
   * @param {Ext.data.Store} store
   */
  onLoad: function(store) {
    this.application.fireEvent('metaboliteload', store);
    this.metabolizable(store.getTotalUnfilteredCount() > 0);
    if (store.getCount() == 1 && !this.getSelectionModel().hasSelection()) {
        Ext.log({}, 'Only one metabolite loaded and its not selected, selecting it');
        this.getSelectionModel().select(0);
    }
    if (store.getTotalUnfilteredCount() === 0) {
        this.showAddStructuresForm();
    }
  },
  /**
   * If metabolite is selected then try to reselect it after load
   * If it fails to reselect it fires a metabolitedeselect event.
   */
  onBeforeLoad: function(store) {
      var me = this;
      var sm = me.getSelectionModel();
      if (sm && sm.hasSelection()) {
          var selected = sm.getSelection()[0].getId();

          var reselect = function() {
              var record = store.getById(selected);
              if (record !== null) {
                sm.select(record);
              } else {
                this.application.fireEvent('metabolitedeselect', selected, 'not found');
              }
              store.removeListener('load', reselect, me);
          };

          store.on('load', reselect , me);
      }
  },
  /**
   * Fetch selection model of metabolite list.
   * @return {Ext.selection.Model}
   */
  getSelectionModel: function() {
    return this.getMetaboliteList().getSelectionModel();
  },
  /**
   * Listens for chromatogram load event.
   * And toggles annotation fieldset in add structures form.
   * @param {Esc.d3.Chromatogram} chromatogram
   */
  onChromatrogramLoad: function(chromatogram) {
    this.hasMSData = chromatogram.data.length > 0;
    this.getMetaboliteAddForm().setDisabledAnnotateFieldset(!this.hasMSData);
  },
  /**
   * Only allow metabolite with a scans to be selected.
   * The extracted ion chromatogram of a metabolite without scans can not be shown because it can not be selected.
   */
  beforeSelect: function(rm, metabolite) {
    return (metabolite.data.nhits > 0);
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
    Ext.log({}, 'Clear metabolite filters');
    this.getMetaboliteList().clearFilters();
    this.application.fireEvent('metabolitenoselect');
  },
  /**
   * Apply scan filter to metabolite store.
   * And shows fragment score column.
   *
   * @param {Number} scanid Scan identifier to filter on.
   */
  applyScanFilter: function(scanid) {
      var store = this.getMetabolitesStore();
      store.sorters.insert(0, new Ext.util.Sorter({
          property: 'score',
          direction: 'DESC'
      }));
      store.setScanFilter(scanid);
      this.getMetaboliteList().showFragmentScoreColumn();
      // TODO in spectra add markers for metabolites present in scan
  },
  /**
   * Removes scan filter from metabolite store.
   * And hides fragment score/deltappm column.
   * And deactivates filters on fragment score/deltappm column if any
   * And resets sort if store is sorted on fragment score/deltappm.
   */
  clearScanFilter: function() {
      var store = this.getMetabolitesStore();
      if ('filters' in store) {
          store.filters.removeAtKey('score');
          store.filters.removeAtKey('deltappm');
      }
      if ('sorters' in store) {
          store.sorters.removeAtKey('score');
          store.sorters.removeAtKey('deltappm');
      }
      store.removeScanFilter();
      this.getMetaboliteList().hideFragmentScoreColumn();
  },
  onPageSizeChange: function(combo) {
      this.getMetabolitesStore().setPageSize(combo.getValue());
  },
  /**
   * Open a new window with metabolites as comma seperated file or sdf.
   * Uses store/proxy/gridfilter state to construct queryString so what you see in grid is what in csv file.
   *
   * @params {String} format Can be csv or sdf.
   */
  download: function(format) {
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
    Ext.apply(params, this.getMetaboliteList().getFilterQuery());

    // visible ordered columns
    params['cols'] = Ext.JSON.encode(this.getMetaboliteList().getVisiblColumnIndices());

    var url = Ext.urlAppend(
        this.application.metabolitesUrl(format),
        Ext.Object.toQueryString(params)
    );
    window.open(url, 'metabolites'+format);
  },
  /**
   * Shortcut do download csv file
   */
  download_csv: function() {
    this.download('csv');
  },
  /**
   * Shortcut do download sdf file
   */
  download_sdf: function() {
    this.download('sdf');
  },
  /**
   * Loads defaults of add structures form and
   * Shows add structures form in metabolites panel.
   */
  showAddStructuresForm: function() {
      this.getMetaboliteAddForm().loadDefaults(this.application.runInfoUrl());
      this.getMetabolitePanel().setActiveItem(1);
  },
  /**
   * Shows list or grid in metabolite panel.
   */
  showGrid: function() {
      this.getMetabolitePanel().setActiveItem(0);
  },
  /**
   * Handler for submit button in Add structures form.
   *
   */
  addStructuresHandler: function() {
    var me = this;
    var form = this.getMetaboliteAddForm().getForm();
    if (form.isValid()) {
      form.submit({
        url: me.application.rpcUrl('add_structures'),
        waitMsg: 'Submitting action ...',
        submitEmptyText: false,
        success: function(fp, o) {
          var response = Ext.JSON.decode(o.response.responseText);
          me.application.fireEvent('rpcsubmitsuccess', response.jobid);
          // switch back to grid so results can still be used while job is running
          me.showGrid();
        },
        failure: function(form, action) {
          if (action.failureType === "server") {
            Ext.Error.raise(Ext.JSON.decode(action.response.responseText));
          } else {
            Ext.Error.raise(action.response.responseText);
          }
        }
      });
    }
  },
  /**
   * Shows metabolize form in modal window
   */
  showMetabolizeForm: function() {
    var me = this;
    if (!this.metabolizeForm) {
        this.metabolizeForm = Ext.create('Esc.magmaweb.view.metabolite.MetabolizeForm');
        this.metabolizeForm.loadDefaults(me.application.runInfoUrl());
    }
    this.metabolizeForm.setDisabledAnnotateFieldset(!this.hasMSData);
    this.metabolizeForm.show();
  },
  metabolizeHandler: function(button) {
    var me = this;
    var wf = this.metabolizeForm;
    var form = wf.getForm();
    if (form.isValid()) {
      form.submit({
        url : me.application.rpcUrl('metabolize'),
        waitMsg: 'Submitting action ...',
        submitEmptyText: false,
        success: function(fp, o) {
          var response = Ext.JSON.decode(o.response.responseText);
          me.application.fireEvent('rpcsubmitsuccess', response.jobid);
          wf.hide();
        },
        failure: function(form, action) {
          wf.hide();
          if (action.failureType === "server") {
            Ext.Error.raise(Ext.JSON.decode(action.response.responseText));
          } else {
            Ext.Error.raise(action.response.responseText);
          }
        }
      });
    }
  },
  /**
   * Shows metabolize form in modal window for one metabolite/structure
   * @param {Ext.data.Model} rec Record to metabolize
   */
  showMetabolizeStructureForm: function(rec) {
    var me = this;
    if (!this.metabolizeStructureForm) {
        this.metabolizeStructureForm = Ext.create('Esc.magmaweb.view.metabolite.MetabolizeOneForm');
        this.metabolizeStructureForm.loadDefaults(me.application.runInfoUrl());
    }
    this.metabolizeStructureForm.setMetabolite(rec);
    this.metabolizeStructureForm.setDisabledAnnotateFieldset(!this.hasMSData);
    this.metabolizeStructureForm.show();
  },
  /**
   * Handler for submit button for metabolize one structure form.
   */
  metabolizeOneHandler: function() {
    var me = this;
    var wf = this.metabolizeStructureForm;
    var form = wf.getForm();
    if (form.isValid()) {
      form.submit({
        url : me.application.rpcUrl('metabolize_one'),
        waitMsg: 'Submitting action ...',
        submitEmptyText: false,
        success: function(fp, o) {
          var response = Ext.JSON.decode(o.response.responseText);
          me.application.fireEvent('rpcsubmitsuccess', response.jobid);
          wf.hide();
        },
        failure: function(form, action) {
          wf.hide();
          if (action.failureType === "server") {
            Ext.Error.raise(Ext.JSON.decode(action.response.responseText));
          } else {
            Ext.Error.raise(action.response.responseText);
          }
        }
      });
    }
  },
  /**
   * Enable or disable metabolize action.
   *
   * @param {Boolean} enabled True is metabolizable
   */
  metabolizable: function(enabled) {
     Ext.getCmp('metabolizeaction').setDisabled(!enabled);
  },
  /**
   * Show actions menu at event xy
   * @param {Ext.Element} tool
   * @param {Ext.EventObject} event
   */
  showActionsMenu: function(tool, event) {
     this.actionsMenu.showAt(event.getXY());
  },
  /**
   * Show download menu at event xy
   * @param {Ext.Element} tool
   * @param {Ext.EventObject} event
   */
  showDownloadMenu: function(tool, event) {
    this.downloadMenu.showAt(event.getXY());
  },
  /**
   * Apply role to user interface.
   * Checks canRun and if false removes all action buttons.
   */
  applyRole: function() {
	  if (this.application.canRun) {
		  return;
	  }
	  this.actionsMenu.getComponent('addstructuresaction').hide();
	  this.actionsMenu.getComponent('metabolizeaction').hide();
	  this.getMetaboliteList().hideCommandsColumn();
  }
});

