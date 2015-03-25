/**
 * Molecule controller.
 *
 * Handles actions performed in molecules views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.controller.Molecules', {
  extend: 'Ext.app.Controller',
  views: [ 'molecule.Panel' ],
  stores: [ 'Molecules' ],
  models: [ 'Molecule' ],
  uses: [
    'Esc.magmaweb.view.molecule.MetabolizeForm',
    'Esc.magmaweb.view.molecule.MetabolizeOneForm'
  ],
  refs: [{
    ref: 'moleculeList', selector: 'moleculelist'
  }, {
    ref: 'moleculeAddForm', selector: 'moleculeaddform'
  }, {
    ref: 'moleculePanel', selector: 'moleculepanel'
  }],
  init: function() {
    Ext.log({}, 'Molecules controller init');
    var me = this;

    // configure store
    var store = this.getMoleculesStore();
    store.pageSize = this.application.getPageSize();
    store.setUrl(this.application.moleculesUrl('json'));
    store.on('load', this.onLoad, this);
    store.on('beforeload', this.onBeforeLoad, this);

    this.control({
      'moleculelist': {
        select: this.onSelect,
        deselect: this.onDeselect,
        beforeselect: this.beforeSelect,
        metabolize: this.showMetabolizeStructureForm
      },
      'moleculelist component[action=pagesizeCombo]': {
        select: this.onPageSizeChange
      },
      'moleculepanel component[action=download]': {
        click: this.showDownloadMenu
      },
      'moleculepanel component[action=actions]': {
        click: this.showActionsMenu
      },
      'moleculepanel component[action=help]': {
          click: this.showHelp
      },
      'moleculeaddform component[action=addstructures]': {
        click: this.addStructuresHandler
      },
      'moleculeaddform component[action=addstructurescancel]': {
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
    this.application.on('peakselect', this.applyMzFilter, this);
    this.application.on('peakdeselect', this.clearMzFilter, this);

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
      me.getMoleculeList().getCommandsColumn().disableAction();
    });

    this.application.addEvents(
        /**
         * @event
         * Triggered when molecule store is loaded.
         * @param {Ext.data.Store} store
         */
        'moleculeload',
        /**
         * @event
         * Triggered when molecule is selected.
         * @param {Number} molid Molecule identifier
         * @param {Esc.magmaweb.model.Molecule} molecule
         */
        'moleculeselect',
        /**
         * @event
         * Triggered when molecule is deselected.
         * @param {Number} molid Molecule identifier
         * @param {Esc.magmaweb.model.Molecule} molecule
         */
        'moleculedeselect',
        /**
         * @event
         * Triggered when molecule selection is cleared
         */
        'moleculenoselect'
    );

    this.actionsMenu = Ext.create('Ext.menu.Menu', {
        items: [{
            iconCls: 'icon-add',
            id: 'addstructuresaction',
            text: 'Add structures',
            listeners: {
                click: {
                    fn: this.showAddStructuresForm,
                    scope: this
                }
            }
        }, {
            text: 'Metabolize',
            id: 'metabolizeaction',
            tooltip: 'Metabolize all structures',
            disabled: true,
            listeners: {
                click: {
                    fn: this.showMetabolizeForm,
                    scope: this
                }
            }
        }, {
            text: 'Clear filters',
            listeners: {
                click: {
                    fn: this.clearFilters,
                    scope: this
                }
            }
        }]
    });

    this.downloadMenu = Ext.create('Ext.menu.Menu', {
        items: [{
            text: 'CSV',
            tooltip: 'Save molecules as comma seperated file',
            listeners: {
                click: {
                    fn: this.download_csv,
                    scope: this
                }
            }
        }, {
            text: 'SDF',
            tooltip: 'Save molecules as sdf',
            listeners: {
                click: {
                    fn: this.download_sdf,
                    scope: this
                }
            }
        }]
    });

    this.application.on('assignmentchanged', function(isAssigned, params) {
        me.getMoleculesStore().load();
    });
  },
  /**
   * Loads molecule store
   */
  onLaunch: function() {
      // store not loaded in init because moleculeload event is fired before listeners of views are registerd
      // the nhits column has an active filter
      // so do not use list.store.load() , but trigger a filter update to load
      this.getMoleculeList().filters.createFilters();
      this.getMoleculesStore().load();
      this.getMoleculeList().setPageSize(this.getMoleculesStore().pageSize);
      this.applyRole();
  },
  /**
   * Listens for molecule store load event.
   * Selects molecule if store only contains 1 molecule.
   *
   * @param {Ext.data.Store} store
   */
  onLoad: function(store) {
    this.application.fireEvent('moleculeload', store);
    this.metabolizable(store.getTotalUnfilteredCount() > 0);
    if (store.getCount() == 1 && !this.getSelectionModel().hasSelection()) {
        Ext.log({}, 'Only one molecule loaded and its not selected, selecting it');
        this.getSelectionModel().select(0);
    }
    if (store.getTotalUnfilteredCount() === 0) {
        this.showAddStructuresForm();
    }
  },
  /**
   * If molecule is selected then try to reselect it after load
   * If it fails to reselect it fires a moleculedeselect event.
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
                this.application.fireEvent('moleculedeselect', selected, 'not found');
              }
              store.removeListener('load', reselect, me);
          };

          store.on('load', reselect , me);
      }
  },
  /**
   * Fetch selection model of molecule list.
   * @return {Ext.selection.Model}
   */
  getSelectionModel: function() {
    return this.getMoleculeList().getSelectionModel();
  },
  /**
   * Listens for chromatogram load event.
   * And toggles annotation fieldset in add structures form.
   * @param {Esc.d3.Chromatogram} chromatogram
   */
  onChromatrogramLoad: function(chromatogram) {
    this.hasMSData = chromatogram.data.length > 0;
    this.getMoleculeAddForm().setDisabledAnnotateFieldset(!this.hasMSData);
  },
  /**
   * Only allow molecule with a scans to be selected.
   * The extracted ion chromatogram of a molecule without scans can not be shown because it can not be selected.
   */
  beforeSelect: function(rm, molecule) {
    return (molecule.data.nhits > 0);
  },
  onSelect: function(rm, molecule) {
    var molid = molecule.data.molid;
    this.application.fireEvent('moleculeselect', molid, molecule);
  },
  onDeselect: function(rm, molecule) {
    var molid = molecule.data.molid;
    this.application.fireEvent('moleculedeselect', molid, molecule);
  },
  /**
   * Remove filters and clears selection
   */
  clearFilters: function() {
    Ext.log({}, 'Clear molecule filters');
    this.getMoleculeList().clearFilters();
    this.application.fireEvent('moleculenoselect');
  },
  /**
   * Apply scan filter to molecule store.
   * And shows fragment score column.
   *
   * @param {Number} scanid Scan identifier to filter on.
   */
  applyScanFilter: function(scanid) {
      var store = this.getMoleculesStore();
      store.sorters.clear();
      store.sorters.addAll([
          new Ext.util.Sorter({
              property: 'score',
              direction: 'ASC'
          }),
          new Ext.util.Sorter({
              property: 'refscore',
              direction: 'DESC'
          }),
          new Ext.util.Sorter({
              property: 'molid',
              direction: 'ASC'
          })
      ]);
      store.setScanFilter(scanid);
      this.getMoleculeList().showFragmentScoreColumn();
  },
  /**
   * Removes scan filter from molecule store.
   * And hides fragment score/deltappm column.
   * And deactivates filters on fragment score/deltappm column if any
   * And resets sort if store is sorted on fragment score/deltappm.
   */
  clearScanFilter: function() {
      var store = this.getMoleculesStore();
      if ('filters' in store) {
          store.filters.removeAtKey('score');
          store.filters.removeAtKey('deltappm');
      }
      if ('sorters' in store) {
          store.sorters.removeAtKey('score');
          store.sorters.removeAtKey('deltappm');
      }
      store.removeScanFilter();
      this.getMoleculeList().hideFragmentScoreColumn();
  },
  applyMzFilter: function(mz, mslevel) {
	  if (mslevel > 1) {
		  return;
	  }
      var store = this.getMetabolitesStore();
	  var list = this.getMetaboliteList();
	  if (store.isFilteredOnScan()) {
		  list.setMzFilterToEqual(mz);
	  }
  },
  clearMzFilter: function(mz, mslevel) {
	  if (mslevel > 1) {
		  return;
	  }
	  var list = this.getMetaboliteList();
	  list.clearMzFilter();
  },
  onPageSizeChange: function(combo) {
      this.getMoleculesStore().setPageSize(combo.getValue());
  },
  /**
   * Open a new window with molecules as comma seperated file or sdf.
   * Uses store/proxy/gridfilter state to construct queryString so what you see in grid is what in csv file.
   *
   * @params {String} format Can be csv or sdf.
   */
  download: function(format) {
    // download needs to make an url with has the same query parameters as the store.load()
    // for load() the store builds an operation object, we need to build this aswell
    var store = this.getMoleculesStore();
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
    Ext.apply(params, this.getMoleculeList().getFilterQuery());

    // visible ordered columns
    params['cols'] = Ext.JSON.encode(this.getMoleculeList().getVisiblColumnIndices());

    var url = Ext.urlAppend(
        this.application.moleculesUrl(format),
        Ext.Object.toQueryString(params)
    );
    window.open(url, 'molecules'+format);
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
   * Shows add structures form in molecules panel.
   */
  showAddStructuresForm: function() {
      this.getMoleculeAddForm().loadDefaults(this.application.runInfoUrl());
      this.getMoleculePanel().setActiveItem(1);
  },
  /**
   * Shows list or grid in molecule panel.
   */
  showGrid: function() {
      this.getMoleculePanel().setActiveItem(0);
  },
  /**
   * Handler for submit button in Add structures form.
   *
   */
  addStructuresHandler: function() {
    var me = this;
    var form = this.getMoleculeAddForm().getForm();
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
        this.metabolizeForm = Ext.create('Esc.magmaweb.view.molecule.MetabolizeForm');
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
   * Shows metabolize form in modal window for one molecule/structure
   * @param {Ext.data.Model} rec Record to metabolize
   */
  showMetabolizeStructureForm: function(rec) {
    var me = this;
    if (!this.metabolizeStructureForm) {
        this.metabolizeStructureForm = Ext.create('Esc.magmaweb.view.molecule.MetabolizeOneForm');
        this.metabolizeStructureForm.loadDefaults(me.application.runInfoUrl());
    }
    this.metabolizeStructureForm.setMolecule(rec);
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
   * Checks run feature and if false removes all action buttons.
   */
  applyRole: function() {
      if (!this.application.features.run) {
          this.actionsMenu.getComponent('addstructuresaction').hide();
          this.actionsMenu.getComponent('metabolizeaction').hide();
          this.getMoleculeList().hideCommandsColumn();
      }
  },
  showHelp: function() {
      this.application.showHelp('molpanel');
  }
});

