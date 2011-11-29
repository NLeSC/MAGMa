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
  setUrl: function(url) {
    this.getProxy().url = url;
  },
  listeners: {
    load: function(store) {
      this.isLoaded = true;
    }
  },
  initComponent: function() {
    console.log('Init Metabolites store');
  },
  /**
   * Removes scan filter from metabolite store
   */
  removeScanFilter: function() {
    // see if already filtered on scanid then remove old filter
    if ('scanid' in this.getProxy().extraParams) {
      delete(this.getProxy().extraParams.scanid);
      this.loadPage(1);
    }
  },
  setScanFilter: function(scanid) {
    this.getProxy().extraParams.scanid = scanid;
    this.loadPage(1);
  }
});

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
    var molcol = Ext.create('Ext.esc.ChemDoodleColumn', {
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
        {text: 'Scans', dataIndex: 'nr_scans', filter: { type: 'numeric' }},
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
  clearFilter: function() {
    this.getView().getFeature('mfilter').clearFilters();
  }
});

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
  },
  onLaunch: function() {
      this.getMetabolitesStore().load();
  },
  onLoad: function(store) {
    this.application.fireEvent('metaboliteload', store);
    if (store.getCount() == 1 && !this.getMetaboliteList().getSelectionModel().hasSelection()) {
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
    console.log('Select metabolite '+metid);
    this.application.fireEvent('metaboliteselect', metid, metabolite);
  },
  onDeselect: function(rm, metabolite) {
    var metid = metabolite.data.metid;
    console.log('Deselect metabolite '+metid);
    this.application.fireEvent('metabolitedeselect', metid, metabolite);
  },
  /**
   * Remove filters and clears selection
   */
  clearFilters: function() {
    console.log('Clear metabolite filter');
    this.getMetaboliteList().clearFilter();
    this.getMetabolitesStore().filter();
    this.application.fireEvent('metabolitenoselect');
  },
  /**
   * If metabolite is selected then try to reselect it after load
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
  applyScanFilter: function(scanid) {
      this.reselectAfterLoad();
      this.getMetabolitesStore().setScanFilter(scanid);
      // TODO in spectra add markers for metabolites present in scan
  },
  clearScanFilter: function() {
      this.reselectAfterLoad();
      this.getMetabolitesStore().removeScanFilter();
  }
});

/**
 * Fragments are loaded when a scan and metabolite are selected.
 */
Ext.define('Esc.msygma.store.Fragments', {
  extend: 'Ext.data.TreeStore',
  model: 'Esc.msygma.model.Fragment',
  autoLoad: false,
  root: { children : [] }, // prevent tree from autoloading
  // TreeStore and Store have different function to fetch record by id, add getById to TreeStore
  getById: function(id) {
    return this.getNodeById(id);
  },
  getNodeByMzMslevel: function(mz, mslevel) {
    return this.getRootNode().findChildBy(function(n) {
      return (n.data.mslevel == mslevel && n.data.mz == mz);
    }, false, true);
  }
});

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
    var fmolcol = Ext.create('Ext.esc.ChemDoodleColumn', {
      pluginId: 'fmolcol',
      text: 'Molecule', dataIndex: 'mol', atomIndex:'atoms',
      canvasClass: 'x-chemdoodle-cols2',
      width: 162,
      initCanvas:function(id, width, height, value,record) {
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
  // Forces all molecules on fragment tree to be drawn
  initMolecules: function() {
    // force fragment molecule rendering, hopyfully canvas have been rendered after spectra has been loaded
    return this.getPlugin('fmolcol').initCanvases();
  }
});

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
    this.application.on('mspectraload', function() { this.getFragmentTree().initMolecules(); }, this);
    this.application.on('peakdeselect', this.clearFragmentSelection, this);
    this.application.on('peakselect', function(mz, mslevel) {
        /**
         * When user selects peak in spectra then select the fragment belonging to peak in fragment tree
         */
        // find fragment based on mz + mslevel
        var node = this.getFragmentsStore().getNodeByMzMslevel(mz, mslevel);
        this.getFragmentTree().getSelectionModel().select([node]);
        if (!node.isLeaf()) {
          if (!node.isExpanded()) {
            node.expand();
          } else {
            this.application.fireEvent('fragmentexpand', node);
          }
        } else {
          // clear product scans
        }
    }, this);
  },
  loadFragments: function (scanid, metid) {
    this.clearFragments();
    console.log('Show fragments of scan '+scanid+' metabolite '+metid);
    var store = this.getFragmentsStore();
    store.setProxy(this.fragmentProxyFactory(scanid, metid));
    store.load();
  },
  // need to change url of fragment proxy so use a factory to create a new proxy foreach scan/metabolite combo
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
  clearFragments: function() {
    console.log('Clearing fragments and mspectra >lvl1');
    this.getFragmentsStore().getRootNode().removeAll();
  },
  onFragmentCollapse: function(fragment) {
    this.application.fireEvent('fragmentcollapse', fragment);
  },
  onFragmentExpand: function(fragment) {
    if (fragment.firstChild == null) {
      return;
    }
    this.application.fireEvent('fragmentexpand', fragment);
  },
  onSelect: function(rm, r) {
    console.log('Selected fragment '+r.id);
    // show child mspectra of selected node or mz
    if (!r.isLeaf()) {
      // onselect then expand
      if (!r.isExpanded()) {
        r.expand();
      } else {
        this.onFragmentExpand(r);
      }
    }
    this.application.fireEvent('fragmentselect', r);
  },
  onDeselect: function(rm, fragment) {
      this.application.fireEvent('fragmentdeselect', fragment);
  },
  clearFragmentSelection: function() {
      this.getFragmentTree().getSelectionModel().deselectAll();
  },
  onLoad: function(t, parent, children) {
    this.application.fireEvent('fragmentload', parent, children);
  }
});

Ext.define('Esc.msygma.resultsApp', {
  extend:'Ext.app.Application',
  constructor: function(config) {
    console.log('Construct app');
    this.initConfig(config);
    this.callParent(arguments);
    return this;
  },
  name: 'Esc.msygma',
  controllers: [ 'Metabolites', 'Fragments', 'Scans' ],
  config: {
    pageSize: 10,
    maxmslevel: 2,
    ms_intensity_cutoff: null,
    urls: {
        fragments: null,
        mspectra: null,
        extractedionchromatogram: null,
        metabolites: null,
        chromatogram: null,
    }
  },
  clearMSpectra: function(mslevel) {
    this.mspectras[mslevel].setData([]);
    this.mspectras[mslevel].scanid = -1;
    Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Scan ... (Level '+mslevel+')');
  },
  loadMSpectra1: function(scanid) {
    this.loadMSpectra(1, scanid, []);
  },
  loadMSpectra2: function(scanid, markers) {
    this.loadMSpectra(2, scanid, markers);
  },
  loadMSpectra: function(mslevel, scanid, markers) {
    var me = this;
    console.log('Loading msspectra level '+mslevel+' with id '+scanid);
    this.mspectras[mslevel].setLoading(true);
    d3.json(Ext.String.format(this.getUrls().mspectra, scanid, mslevel), function(data) {
      if (!data) {
        Ext.MessageBox.show({
          title: 'Unable to find scan',
          msg: 'Level '+mslevel+' scan with id '+scanid+' was not found',
          buttons: Ext.MessageBox.OK,
          icon: Ext.MessageBox.ERROR
        });
        me.mspectras[mslevel].setLoading(false);
        return;
      }
      Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Scan '+scanid+' (Level '+mslevel+')');
      me.mspectras[mslevel].setLoading(false);
      me.mspectras[mslevel].scanid = scanid;
      me.mspectras[mslevel].cutoff = data.cutoff;
      me.mspectras[mslevel].setData(data.peaks);
      me.mspectras[mslevel].setMarkers(markers);

      me.fireEvent('mspectraload', scanid, mslevel);
    });
  },
  loadMSpectraFromFragment: function(fragment) {
      // on expand load child mspectra if needed
      if (fragment.firstChild.data.scanid != this.mspectras[fragment.firstChild.data.mslevel].scanid) {
        this.loadMSpectra(
          fragment.firstChild.data.mslevel,
          fragment.firstChild.data.scanid,
          fragment.childNodes.map(function(r) { return {mz: r.data.mz}; })
        );
      }
  },
  selectPeak: function(fragment) {
      // select peak belonging to r
      this.mspectras[fragment.data.mslevel].selectPeak(fragment.data.mz);
      // select peaks of parents of fragment in parent scans
      if (fragment.data.mslevel==1) {
        this.mspectras[1].selectPeak(fragment.data.mz);
        for (var i = 2; i <= this.getMaxmslevel(); i++) {
          this.mspectras[i].clearPeakSelection();
        }
      } else if (fragment.data.mslevel==2) {
        for (var i = 3; i <= this.getMaxmslevel(); i++) {
          this.mspectras[i].clearPeakSelection();
        }
        this.mspectras[2].selectPeak(fragment.data.mz);
        this.mspectras[1].selectPeak(fragment.parentNode.data.mz);
      } else if (fragment.data.mslevel>=3) {
        // TODO make selecting parent Node work for mslevel>3
        this.mspectras[2].selectPeak(fragment.parentNode.data.mz);
        this.mspectras[1].selectPeak(fragment.parentNode.parentNode.data.mz);
      }
  },
  deselectPeak: function(fragment) {
      this.mspectras[fragment.data.mslevel].clearPeakSelection();
      // also unselect peaks in child scans
      for (var i = fragment.data.mslevel+1; i <= this.getMaxmslevel(); i++) {
          this.mspectras[i].clearPeakSelection();
      }
  },
  /**
   * Depending on which mslevel the fragment was found and if has children
   * MSpectra are loaded and peaks selected.
   */
  loadMSpectras: function(parent, children) {
      // show peaks in lvl1 scan
      if ('id' in parent.data && parent.data.id == 'root') {
        console.log('Loaded metabolite as fragment');
        // add mz of metabolites as markers to lvl1 scan
        this.mspectras[1].setMarkers(
          children.map(function(r) { return {mz: r.data.mz}; })
        );
        var metabolite_fragment = children[0];
        this.mspectras[1].selectPeak(metabolite_fragment.data.mz);
        if (metabolite_fragment.hasChildNodes()) {
          this.loadMSpectra2(
            metabolite_fragment.childNodes[0].data.scanid,
            metabolite_fragment.childNodes.map(function(r) { return {mz: r.data.mz}; })
          );
        }
      } else if (parent.data.mslevel == 1) {
        console.log('Loaded lvl2 fragments of metabolite ');
        // load the scan of first child
        // add mz of metabolites as markers to lvl2 scan
        this.loadMSpectra2(
          children[0].data.scanid,
          children.map(function(r) { return {mz: r.data.mz}; })
        );
        this.mspectras[1].selectPeak(parent.data.mz);
      } else if (parent.data.mslevel >= 2) {
        console.log('Loaded lvl'+(parent.data.mslevel+1)+' fragments of metabolite');
        // select parent peak in lvl2 mspectra
        this.mspectras[2].selectPeak(parent.data.mz);
        // TODO select parent peaks if n.data.mslevel>2
      }
  },
  launch: function() {

    // uncomment to see all application events fired in console
    Ext.util.Observable.capture(this, function() { console.log(arguments);return true;});

    this.on('selectscan', this.loadMSpectra1, this);
    // when a metabolite and scan are selected then load fragments
    this.selected = { scanid: false, metid: false };
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
        me.clearMSpectra(1);
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
    this.on('scanandmetabolitenoselect', function() {
        this.mspectras[1].setMarkers([]);
        for (var i = 2; i <= this.maxmslevel; i++) {
          this.clearMSpectra(i);
        }
    }, this);
    this.on('fragmentcollapse', function(fragment) {
        this.mspectras.forEach(function(ms, i) {
            if (i > fragment.data.mslevel ) {
              this.clearMSpectra(i);
            }
        }, this);
    }, this);
    this.on('fragmentexpand', this.loadMSpectraFromFragment, this);
    this.on('scanandmetabolitenoselect', function() {
        this.mspectras[1].setMarkers([]);
        for (var i = 2; i <= this.maxmslevel; i++) {
          this.clearMSpectra(i);
        }
    }, this);
    /**
     * When user selects fragment in tree then select the peak in the mspectra
     */
    this.on('fragmentselect', this.selectPeak, this);
    this.on('fragmentdeselect', this.deselectPeak, this);
    this.on('fragmentload', this.loadMSpectras, this);

    console.log('Launch app');
    var me = this;
    var config = me.config;

    var msspectrapanels = [];
    this.mspectras = [];
    for (var mslevel = 1; mslevel <= this.getMaxmslevel(); mslevel++) {
      this.mspectras[mslevel] = Ext.create('Ext.esc.MSpectra', {
        mslevel: mslevel,
        emptyText: (
            mslevel==1 ?
            'Select a scan in the chromatogram' :
            'Select a fragment to show its level '+mslevel+' scan'
        ),
        listeners: {
          selectpeak: function(mz) {
            me.fireEvent('peakselect', mz, this.mslevel);
            // TODO unselect peaks of child scans
          },
          unselectpeak: function(mz) {
            me.fireEvent('peakdeselect', mz, this.mslevel);
          }
        }
      });
      msspectrapanels.push({
        title: 'Scan ... (Level '+mslevel+')',
        id: 'mspectra'+mslevel+'panel',
        collapsible: true,
        items: this.mspectras[mslevel]
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
