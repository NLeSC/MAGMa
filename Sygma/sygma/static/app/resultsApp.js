/**
 * Metabolite model.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.msygma.model.Metabolite', {
  idProperty: 'metid',
  extend:'Ext.data.Model',
  fields: [{
    name: 'metid',
  },{
    name: 'mol', type:'string'
  },{
    name: 'level'
  },{
    name: 'probability'
  },{
    name: 'reactionsequence'
  },{
    name: 'smiles'
  },{
    name: 'molformula'
  },{
    name: 'isquery', type: 'bool'
  },{
    name: 'origin'
  },{
    name: 'nhits'
  },{
    name: 'atoms', defaultValue: [] // array of atom indexes of molecule which are the substructure of the query
  },{
    name: 'nr_scans', type:'number'
  },{
    name: 'scans', defaultValue: [] // Filled when metabolite is selected
  }]
});

/**
 * Fragment model.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.msygma.model.Fragment', {
  extend:'Ext.data.Model',
  idProperty: 'fragid',
  fields: [{
    name: 'fragid',
  },{
    name: 'metid',
  },{
    name: 'scanid'
  },{
    name: 'mz'
  },{
    name: 'mass'
  },{
    name: 'score'
  },{
    name: 'parentfragid'
  },{
    name: 'mol', type: 'string'
  },{
    name: 'atoms', defaultValue: [] // array of atom indexes of molecule which are the substructure of the query
  },{
    name: 'deltah'
  },{
    name: 'mslevel'
  }],
  hasMany: { model: 'Fragment', name:'children' }
});

/**
 * Metabolite store.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.msygma.store.Metabolites', {
  extend: 'Ext.data.Store',
  model: 'Esc.msygma.model.Metabolite',
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

/**
 * Grid of metabolites.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.msygma.view.metabolite.List', {
  extend: 'Ext.grid.Panel',
  alias: 'widget.metabolitelist',
  title: 'Query molecules & Metabolites',
  store: 'Metabolites',
  selModel: Ext.create('Ext.selection.CheckboxModel', {
    allowDeselect: true,
    mode: 'SINGLE'
  }),
  scroll: false,
  viewConfig: {
    autoScroll: true,
  },
  dockedItems: [{
    xtype: 'pagingtoolbar',
    store: 'Metabolites',   // same store GridPanel is using
    dock: 'bottom',
    displayInfo: true,
    items: [{
      text: 'Clear filters',
      action: 'clear'
    }]
  }],
  initComponent: function() {
    console.log('Init met grid');
    var molcol = Ext.create('Esc.chemdoodle.Column', {
      text: 'Molecule', dataIndex: 'mol',
      width: 162
    });

    var mfilters = Ext.create('Ext.ux.grid.FiltersFeature',{
      id: 'mfilter',
      encode: true,
    });

    Ext.apply(this, {
      columns: [
        {text: 'ID', dataIndex: 'metid', hidden: true},
        molcol,
        {text: 'Level', dataIndex: 'level', filter: { type: 'list',  options: ['0','1','2','3'] }, hidden:true},
        {text: 'Probability', dataIndex: 'probability', filter: { type: 'numeric' }},
        {text: 'Reaction seq.', dataIndex: 'reactionsequence', flex:1, filter: { type: 'string' }, renderer: function(v) {
          return '<ol><li>'+v.replace("\n","</li><li>")+'</li></ol>';
        }},
        {text: 'Scans', dataIndex: 'nr_scans', filter: { type: 'numeric'
            , value:{gt:0}, active: true
        }},
        {text: 'Smile', dataIndex: 'smiles', hidden:true},
        {text: 'Formula', dataIndex: 'molformula', filter: { type: 'string' }},
        {text: 'Query', dataIndex: 'isquery', xtype:'booleancolumn', trueText:'Yes', falseText:'No', filter: { type: 'boolean' }},
        {text: 'Name', dataIndex: 'origin', hidden: true, filter: { type: 'string' }},
      ],
      plugins: [molcol],
      features: [mfilters],
    });
    this.callParent(arguments);
  },
  /**
   * Clears all filters applied to metabolites
   */
  clearFilters: function() {
    this.getView().getFeature('mfilter').clearFilters();
  }
});

/**
 * Metabolite controller.
 *
 * Handles actions performed in metabolites views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.msygma.controller.Metabolites', {
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

    var grid = this.getMetaboliteListView();
    grid.pageSize = this.application.getPageSize();

    this.control({
      'metabolitelist': {
        select: this.onSelect,
        deselect: this.onDeselect,
        beforeselect: this.beforeSelect,
      },
      'metabolitelist button[action=clear]': {
        click: this.clearFilters
      },
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
         * @param {Esc.msygma.model.Metabolite} metabolite
         */
        'metaboliteselect',
        /**
         * @event
         * Triggered when metabolite is deselected.
         * @param {Number} metid Metabolite identifier
         * @param {Esc.msygma.model.Metabolite} metabolite
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
      // not loaded in init because metaboliteload event is fired before listeners are registerd
      // the nr_scans column has an active filter
      // so do not use list.store.load() , but trigger a filter update to load
      this.getMetaboliteList().filters.createFilters();
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
   * Tries to keep selection.
   *
   * @param {Number} scanid Scan identifier to filter on.
   */
  applyScanFilter: function(scanid) {
      this.reselectAfterLoad();
      this.getMetabolitesStore().setScanFilter(scanid);
      // TODO in spectra add markers for metabolites present in scan
  },
  /**
   * Removes scan filter from metabolite store.
   */
  clearScanFilter: function() {
      this.getMetabolitesStore().removeScanFilter();
  }
});

/**
 * Store for fragments.
 *
 * Fragments are loaded when a scan and metabolite are selected.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.msygma.store.Fragments', {
  extend: 'Ext.data.TreeStore',
  model: 'Esc.msygma.model.Fragment',
  autoLoad: false,
  root: { children : [] }, // prevent tree from autoloading
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

/**
 * Fragment tree.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.msygma.view.fragment.Tree', {
  extend: 'Ext.tree.Panel',
  alias: 'widget.fragmenttree',
  title: 'Fragments',
  store: 'Fragments',
  selType: 'checkboxmodel',
  cls: 'fragmenttree', // So height of column can be set with .fragmenttree .x-grid-cell-inner{height: 106px !important;}
  multiSelect: false,
  rootVisible: false,
  singleExpand: true,
  scroll: false,
  viewConfig: {
    // animate is default true causing refresh event to be blocked
    // we use refresh event to render molecules
    // so after expanding a node the refresh was not fired causing all prev. rendered mols to disappear
    // now we turn off animate, so refresh events are fired and mols can be rendered
    animate: false,
    autoScroll: true,
    blockRefresh: false,
    emptyText: 'Select a metabolite and scan, to show its fragments',
  },
  initComponent: function() {
    console.log('Init fragment tree');

    // atoms property is array filled with fragment atoms that need to be black
    // bonds having both atoms in array are black
    var hlspecs = new ChemDoodle.structures.VisualSpecifications();
    hlspecs.bonds_color = 'black';
    hlspecs.atoms_color = 'black';
    var fmolcol = Ext.create('Esc.chemdoodle.Column', {
      pluginId: 'fmolcol',
      text: 'Molecule', dataIndex: 'mol', atomIndex:'atoms',
      canvasClass: 'x-chemdoodle-cols2',
      width: 162,
      initCanvas: function(id, width, height, value,record) {
        var c = new ChemDoodle.ViewerCanvas(id, width, height,true);
        c.specs.bonds_color = 'cyan';
        c.specs.atoms_color = 'cyan';
        var m = ChemDoodle.readMOL(value);
        var fragmentAtoms = record.data[this.atomIndex].split(',');
        m.atoms.forEach(function(v, i) {
          if (fragmentAtoms.indexOf(i+'') != -1) {
            // use highlight visual spec for atoms which are part of fragment
            v.specs = hlspecs;
            // also highlight all bonds connected to fragment atoms
            m.getBonds(v).forEach(function(b) {
              // only color bond black if both atoms are in fragmentAtoms
              var ni = m.atoms.indexOf(b.getNeighbor(v));
              if (ni != -1 && fragmentAtoms.indexOf(ni+'') != -1) {
                b.specs = hlspecs;
              }
            });
          }
        });
        c.loadMolecule(m);
      }
    });

    Ext.apply(this, {
      columns: [
        { text: 'Score', dataIndex: 'score', xtype: 'treecolumn', width: 120},
        fmolcol,
        { text: 'ID', dataIndex: 'fragid', hidden: true},
        { text: 'Scan', dataIndex: 'scanid', hidden: false},
        { text: 'Metabolite', dataIndex: 'metid', hidden: true},
        { text: 'M/z', dataIndex: 'mz'},
        { text: 'Mass', dataIndex: 'mass'},
        { text: 'Level', dataIndex: 'mslevel'},
        { text: 'Fragment atoms', dataIndex: 'atoms', hidden: true},
        { text: 'H Delta', dataIndex: 'deltah'}
      ],
      plugins: [fmolcol],
    });

    this.callParent(arguments);
  },
  /**
   * Forces all molecules on fragment tree to be drawn
   */
  initMolecules: function() {
    // force fragment molecule rendering, hopyfully canvas have been rendered after spectra has been loaded
    return this.getPlugin('fmolcol').initCanvases();
  }
});

/**
 * Fragment controller.
 *
 * Handles actions performed on the fragment views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.msygma.controller.Fragments', {
  extend: 'Ext.app.Controller',
  stores: [ 'Fragments' ],
  models: [ 'Fragment' ],
  views: [ 'fragment.Tree' ],
  refs: [{
    ref: 'fragmentTree', selector: 'fragmenttree'
  }],
  init: function() {
    console.log('Fragments controller init');

    this.getFragmentsStore().on('load', this.onLoad, this);

    this.control({
      'fragmenttree': {
        select: this.onSelect,
        deselect: this.onDeselect,
        itemcollapse: this.onFragmentCollapse,
        itemexpand: this.onFragmentExpand
      }
    });

    this.application.on('scanandmetaboliteselect', this.loadFragments, this);
    this.application.on('scanandmetabolitenoselect', this.clearFragments, this);
    this.application.on('mspectraload', this.initMolecules, this);
    this.application.on('peakdeselect', this.clearFragmentSelection, this);
    this.application.on('peakselect', this.selectFragmentByPeak, this);

    this.application.addEvents(
      /**
       * @event
       * Triggered when a fragment node is collapsed.
       * @param {Esc.msygma.model.Fragment} fragment Fragment which has been collapsed.
       */
      'fragmentcollapse',
      /**
       * @event
       * Triggered when a fragment node is expanded.
       * @param {Esc.msygma.model.Fragment} fragment Fragment which has been expanded.
       */
      'fragmentexpand',
      /**
       * @event
       * Triggered when a fragment node is selected.
       * @param {Esc.msygma.model.Fragment} fragment Fragment which has been selected.
       */
      'fragmentselect',
      /**
       * @event
       * Triggered when a fragment node is deselected.
       * @param {Esc.msygma.model.Fragment} fragment Fragment which has been deselected.
       */
      'fragmentdeselect',
      /**
       * @event
       * Triggered when a children of a fragment node are loaded.
       * @param {Esc.msygma.model.Fragment} parent
       * @param {Array} children Array of fragment children.
       */
      'fragmentload'
    );
  },
  /**
   * Loads lvl 1 and 2 fragments of a metabolite scan combination.
   *
   * @param {Number} scanid Scan identifier.
   * @param {Number} metid Metabolite idenfitier.
   */
  loadFragments: function (scanid, metid) {
    this.clearFragments();
    console.log('Show fragments of scan '+scanid+' metabolite '+metid);
    var store = this.getFragmentsStore();
    store.setProxy(this.fragmentProxyFactory(scanid, metid));
    store.load();
  },
  /**
   * Need to change url of fragment proxy so use a factory to create a new proxy foreach scan/metabolite combo
   *
   * @param {Number} scanid Scan identifier.
   * @param {Number} metid Metabolite idenfitier.
   * @private
   */
  fragmentProxyFactory: function (scanid, metid) {
    return Ext.create('Ext.data.proxy.Ajax', {
      // url is build when scan and metabolite are selected
      url: Ext.String.format(this.application.getUrls().fragments, scanid, metid),
      reader: {
          type: 'json',
          root: 'children',
          idProperty: 'fragid'
      }
    });
  },
  /**
   * Clears fragments from store.
   */
  clearFragments: function() {
    console.log('Clearing fragments and mspectra >lvl1');
    this.getFragmentsStore().getRootNode().removeAll();
  },
  onFragmentCollapse: function(fragment) {
    this.application.fireEvent('fragmentcollapse', fragment);
  },
  onFragmentExpand: function(fragment) {
    if (fragment.firstChild == null) {
      return; // root node auto expands, but is no fragment, so dont fire event
    }
    this.application.fireEvent('fragmentexpand', fragment);
  },
  onSelect: function(rm, r) {
    console.log('Selected fragment '+r.id);
    // show child mspectra of selected node or mz
    if (!r.isLeaf()) {
      // onselect then expand
      if (r.isExpanded()) {
        this.onFragmentExpand(r);
      } else {
        r.expand();
      }
    }
    this.application.fireEvent('fragmentselect', r);
  },
  onDeselect: function(rm, fragment) {
      this.application.fireEvent('fragmentdeselect', fragment);
  },
  /**
   * Clears fragment selection.
   */
  clearFragmentSelection: function() {
      this.getFragmentTree().getSelectionModel().deselectAll();
  },
  onLoad: function(t, parent, children) {
    this.application.fireEvent('fragmentload', parent, children);
  },
  selectFragment: function(fragment) {
    this.getFragmentTree().getSelectionModel().select([fragment]);
    if (!fragment.isLeaf()) {
      if (fragment.isExpanded()) {
        this.application.fireEvent('fragmentexpand', fragment);
      } else {
        fragment.expand();
      }
    }
  },
  /**
   * When user selects peak in spectra then select the fragment belonging to peak in fragment tree
   *
   * @param {Number} mz m/z of peak
   * @param {Number} mslevel MS level of peak
   */
  selectFragmentByPeak: function(mz, mslevel) {
    // find fragment based on mz + mslevel
    var node = this.getFragmentsStore().getNodeByMzMslevel(mz, mslevel);
    this.selectFragment(node);
  },
  /**
   * Forces molecules canvases to be drawn
   */
  initMolecules: function() {
    this.getFragmentTree().initMolecules();
  }
});

/**
 * MSygma results application
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.msygma.resultsApp', {
  extend:'Ext.app.Application',
  constructor: function(config) {
    console.log('Construct app');
    this.initConfig(config);
    this.callParent(arguments);
    return this;
  },
  name: 'Esc.msygma',
  controllers: [ 'Metabolites', 'Fragments', 'Scans', 'MSpectras' ],
  config: {
    /**
     * Metabolite grid page size.
     * @cfg {Number}
     */
    pageSize: 10,
    /**
     * Maximum MS level or nr of MS levels.
     * @cfg {Number}
     */
    maxmslevel: 2,
    /**
     * MS intensity cutoff. Intensity under which peaks are ignored.
     * @cfg {Number}
     */
    ms_intensity_cutoff: null,
    /**
     * Endpoints/templates for contacting server.
     * @cfg {Object}
     */
    urls: {
        /**
         * Fragments endpoint.
         * Tokenized string with scanid and metid tokens.
         * @cfg {String} urls.fragments
         */
        fragments: null,
        /**
         * MSpectra endpoint.
         * Tokenized string with scanid and mslevel tokens.
         * @cfg {String} urls.mspectra
         */
        mspectra: null,
        /**
         * Extracted ion chromatogram endpoint.
         * Tokenized string with metid token.
         * @cfg {String} urls.extractedionchromatogram
         */
        extractedionchromatogram: null,
        /**
         * Metabolites endpoint.
         * @cfg {String} urls.metabolites
         */
        metabolites: null,
        /**
         * Chromatogram endpoint.
         * @cfg {String} urls.chromatogram
         */
        chromatogram: null
    }
  },
  /**
   * when a metabolite and scan are selected then load fragments
   * @property {Object} selected
   * @property {Boolean} selected.scanid Scan identifier
   * @property {Boolean} selected.metid Metabolite identifier
   */
  selected: { scanid: false, metid: false },
  /**
   * Creates mspectraspanels and viewport and fires/listens for mspectra events
   */
  launch: function() {
    this.addEvents(
      /**
       * @event
       * Triggered when a metabolite and scan are selected together.
       * @param {Number} scanid Scan identifier.
       * @param {Number} metid Metabolite identifier.
       */
      'scanandmetaboliteselect',
      /**
       * @event
       * Triggered when a metabolite and scan are no longer selected together.
       * @param {Number} scanid Scan identifier.
       * @param {Number} metid Metabolite identifier.
       */
      'scanandmetabolitenoselect'
    );

    // uncomment to see all application events fired in console
    Ext.util.Observable.capture(this, function() { console.log(arguments);return true;});

    this.on('metaboliteselect', function(metid) {
      this.selected.metid = metid;
      if (this.selected.metid && this.selected.scanid) {
        this.fireEvent('scanandmetaboliteselect', this.selected.scanid, metid);
      }
    }, this);
    this.on('selectscan', function(scanid) {
        this.selected.scanid = scanid;
        if (this.selected.metid && this.selected.scanid) {
            this.fireEvent('scanandmetaboliteselect', scanid, this.selected.metid);
        }
    }, this);
    this.on('noselectscan', function() {
        if (this.selected.metid && this.selected.scanid) {
            this.fireEvent('scanandmetabolitenoselect');
        }
        this.selected.scanid = false;
    }, this);
    this.on('metabolitedeselect', function() {
        if (this.selected.metid && this.selected.scanid) {
            this.fireEvent('scanandmetabolitenoselect');
        }
        this.selected.metid = false;
    }, this);
    this.on('metabolitenoselect', function() {
        if (this.selected.metid && this.selected.scanid) {
            this.fireEvent('scanandmetabolitenoselect');
        }
        this.selected.metid = false;
    }, this);

    this.on('mspectraload', function(scanid, mslevel) {
      Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Level '+mslevel+' scan '+scanid);
    });
    this.on('mspectraclear', function(mslevel) {
      Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Level '+mslevel+' scan ...');
    });
    this.on('peakmouseover', function(peak, mslevel, scanid) {
      Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Level '+mslevel+' scan '+scanid+' (m/z='+peak.mz+', intensity='+peak.intensity+')');
    });

    console.log('Launch app');

    var msspectrapanels = [];
    var mspectras = this.getController('MSpectras').mspectras;
    for (var mslevel = 1; mslevel <= this.getMaxmslevel(); mslevel++) {
      msspectrapanels.push({
        title: 'Level '+mslevel+' scan ...',
        id: 'mspectra'+mslevel+'panel',
        collapsible: true,
        items: mspectras[mslevel]
      });
    }

    var master_side = Ext.create('Ext.panel.Panel', {
      // master side
      region: 'center',
      layout: 'border',
      border: false,
      items:[{
        region:'center',
        border: false,
        xtype: 'metabolitelist'
      },{
        region:'south',
        hideCollapseTool: true,
        collapsible: true,
        height: '50%',
        split: true,
        xtype: 'scanchromatogram',
        border: false,
      }]
    });

    // detail side
    var detail_side = Ext.create('Ext.panel.Panel', {
      region: 'east',
      split: true,
      collapsible: true,
      layout: 'border',
      width: 600,
      hideCollapseTool: true,
      border: false,
      items:[{
        region: 'center',
        xtype: 'fragmenttree',
        border: false
      },{
        region:'south',
        height: '50%',
        layout: 'border',
        split: true,
        collapsible: true,
        hideCollapseTool: true,
        border: false,
        items:[{
          id: 'mspectrapanel',
          region: 'center',
          layout: {
            type: 'vbox',
            align: 'stretch'
          },
          border: false,
          defaults: {
            flex: 1,
            layout:'fit',
            border: false
          },
          items: msspectrapanels
        }],
      }]
    });

    Ext.create('Ext.Viewport', {
      layout: 'border',
      items:[ master_side, detail_side ]
    });
  }
});
