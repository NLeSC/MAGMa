/**
 * Metabolite controller.
 *
 * Handles actions performed in metabolites views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.controller.Metabolites', {
  extend: 'Ext.app.Controller',
  views: [ 'metabolite.List' ],
  stores: [ 'Metabolites' ],
  models: [ 'Metabolite' ],
  uses: [
    'Ext.window.Window',
    'Ext.form.Panel',
    'Ext.form.field.Hidden',
    'Ext.form.field.Display',
    'Esc.magmaweb.view.metabolite.AddFieldSet',
    'Esc.magmaweb.view.metabolite.MetabolizeFieldSet',
    'Esc.magmaweb.view.fragment.AnnotateFieldSet'
  ],
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
        beforeselect: this.beforeSelect,
        metabolize: this.showMetabolizeStructureForm
      },
      'metabolitelist component[action=pagesize]': {
        select: this.onPageSizeChange
      },
      'metabolitelist component[action=download]': {
        click: this.download
      },
      'metabolitelist component[action=actions]': {
        click: this.showActionsMenu
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
    this.metabolizable(store.getTotalCount() > 0);
  },
  /**
   * Listens for chromatogram load event.
   * @param {Esc.d3.Chromatagram} chromatogram
   */
  onChromatrogramLoad: function(chromatogram) {
    this.hasMSData = chromatogram.data.length > 0;
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
  },
  /**
   * Shows add structures form in modal window.
   */
  showAddStructuresForm: function() {
      var me = this;
      if (!this.addStructuresForm) {
          this.addStructuresForm = Ext.create('Ext.window.Window', {
              title: 'Add structure(s)',
              height: 500,
              width: 600,
              layout: 'fit',
              modal: true,
              closeAction: 'hide',
              items: {
                    xtype: 'form',
                    bodyPadding: 5,
                    defaults: { bodyPadding: 5 },
                    border: false,
                    autoScroll: true,
                    url: me.application.rpcUrl('add_structures'),
                    items : [{
                        xtype : 'addstructurefieldset'
                    }, {
                        xtype : 'metabolizefieldset',
                        checkboxToggle: true,
                        checkboxName: 'metabolize',
                        collapsed : true,
                        collapsible : true
                    }, {
                        xtype : 'annotatefieldset',
                        disabled: !this.hasMSData,
                        collapsed : true,
                        collapsible : true
                    }],
                    buttons: [{
                        text: 'Submit',
                        handler: function() {
                            var form = this.up('form');
                            me.actionHandler(form);
                        }
                    }, {
                        text: 'Reset',
                        handler: function() {
                            this.up('form').getForm().reset();
                        }
                    }]
              }
          });
      }
      this.addStructuresForm.show();
  },
  /**
   * Shows metabolize form in modal window
   */
  showMetabolizeForm: function() {
    var me = this;
    if (!this.MetabolizeForm) {
        this.MetabolizeForm = Ext.create('Ext.window.Window', {
            title: 'Metabolize all structures',
            modal: true,
            height: 300,
            width: 600,
            layout: 'fit',
            closeAction: 'hide',
            items: {
                xtype: 'form',
                bodyPadding: 5,
                defaults: { bodyPadding: 5 },
                border: false,
                autoScroll: true,
                url: me.application.rpcUrl('metabolize'),
                items: [{
                    xtype : 'metabolizefieldset'
                }, {
                    xtype : 'annotatefieldset',
                    disabled: !this.hasMSData,
                    collapsed : true,
                    collapsible : true
                }],
                buttons: [{
                    text: 'Submit',
                    handler: function() {
                        var form = this.up('form');
                        me.actionHandler(form);
                    }
                }, {
                    text: 'Reset',
                    handler: function() {
                        this.up('form').getForm().reset();
                    }
                }]
            }
        });
    }
    this.MetabolizeForm.show();
  },
  /**
   * Shows metabolize form in modal window for one metabolite/structure
   * @param {Ext.data.Model} rec Record to metabolize
   */
  showMetabolizeStructureForm: function(rec) {
    var me = this;
    if (!this.metabolizeStructureForm) {
        this.metabolizeStructureForm = Ext.create('Ext.window.Window', {
            title: 'Metabolize',
            modal: true,
            height: 300,
            width: 600,
            layout: 'fit',
            closeAction: 'hide',
            items: {
                xtype: 'form',
                url: me.application.rpcUrl('metabolize_one'),
                bodyPadding: 5,
                defaults: { bodyPadding: 5 },
                border: false,
                autoScroll: true,
                items: [{
                    xtype: 'displayfield',
                    fieldLabel: 'Name',
                    value: rec.get('origin')
                },{
                    xtype: 'hiddenfield',
                    name: 'metid',
                    value: rec.get('metid')
                },{
                    xtype : 'metabolizefieldset'
                }, {
                    xtype : 'annotatefieldset',
                    disabled: !this.hasMSData,
                    collapsed : true,
                    collapsible : true
                }],
                buttons: [{
                    text: 'Submit',
                    handler: function() {
                        var form = this.up('form');
                        me.actionHandler(form);
                    }
                }, {
                    text: 'Reset',
                    handler: function() {
                        this.up('form').getForm().reset();
                    }
                }]
            }
        });
    }
    this.metabolizeStructureForm.show();
  },
  metabolizable: function(enabled) {
     Ext.getCmp('metabolizeaction').setDisabled(!enabled);
  },
  actionHandler: function(form) {
      var wf = form.up('window');
      form = form.getForm();
      if (form.isValid()) {
          form.submit({
              waitMsg: 'Submitting action ...',
              success: function(fp, o) {
                  console.log('Action submitted');
                  wf.hide();
              },
              failure: function(form, action) {
                  console.log(action.failureType);
                  console.log(action.result);
                  wf.hide();
              }
          });
      }
  },
  showActionsMenu: function(tool, event) {
    if (this.actionsMenu.isHidden()) {
        this.actionsMenu.showAt(event.getXY());
    } else {
        this.actionsMenu.hide();
    }
  }
});

