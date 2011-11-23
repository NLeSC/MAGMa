<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <title>MSygma - Results</title>
  <meta http-equiv="Content-Type" content="text/html;charset=UTF-8"/>
<link rel="stylesheet" href="${request.static_url('sygma:static/ChemDoodleWeb/ChemDoodleWeb.css')}" type="text/css"></link>
<link rel="stylesheet" href="${request.static_url('sygma:static/ext-4.0.7-gpl/resources/css/ext-all.css')}" type="text/css"></link>
<link rel="stylesheet" type="text/css" href="${request.static_url('sygma:static/ext-4.0.7-gpl/examples/ux/grid/css/GridFilters.css')}" />
<link rel="stylesheet" type="text/css" href="${request.static_url('sygma:static/ext-4.0.7-gpl/examples/ux/grid/css/RangeMenu.css')}" />
<script type="text/javascript" src="${request.static_url('sygma:static/ext-4.0.7-gpl/bootstrap.js')}"></script>
<script type="text/javascript" src="${request.static_url('sygma:static/ChemDoodleWeb/ChemDoodleWeb-libs.js')}"></script>
<script type="text/javascript" src="${request.static_url('sygma:static/ChemDoodleWeb/ChemDoodleWeb.js')}"></script>
<script type="text/javascript" src="${request.static_url('sygma:static/d3/d3.js')}"></script>
<style type="text/css">

path.line {
  fill: none;
  stroke: #66d;
  stroke-width: 2px;
}

path.metaboliteline {
  fill: none;
  stroke: #6d6;
  stroke-width: 2px;
}

.axis {
  shape-rendering: crispEdges;
}

.y.axis line, .y.axis path, .x.axis line, .x.axis path {
  fill: none;
  stroke: #ccc;
}

.peaks {
  fill: none;
  stroke-width: 1px;
  stroke: "#111";
}

.cutoffline {
  fill: none;
  stroke-width: 1px;
  stroke: #ddd;
  stroke-dasharray: 5px, 5px;
}

line.peak {
  stroke-width: 2px;
  stroke: #eee;
}

path.marker {
  stroke: lightgreen;
  fill: none;
  stroke-width: 1.5px;
  cursor: pointer;
}

.selectedscan, .selectedpeak {
  stroke: darkgreen !important;
  fill: darkgreen !important;
}

line.mspeak {
  stroke-width: 1px;
  stroke: black;
}

.fragmenttree .x-grid-cell-inner{
    height: 106px !important;
}

</style><script type="text/javascript">

Ext.Loader.setConfig({
  enabled: true,
  paths: {
    'Ext': '${request.static_url('sygma:static/')}',
    'Ext.ux': '${request.static_url('sygma:static/ext-4.0.7-gpl/examples/ux')}'
  }
});

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
  }
});

Ext.define('Esc.msygma.view.metabolite.List', {
  extend: 'Ext.grid.Panel',
  alias: 'widget.metabolitelist',
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

    // configure store
    var store = this.getMetabolitesStore();
    store.pageSize = this.application.getPageSize();
    store.setUrl(this.application.getUrls().metabolites);
    store.on('load', this.onLoad);
    store.load();

    var grid = this.getMetaboliteListView();
    grid.pageSize = this.application.getPageSize();
    console.log('Metabolites controller init');

    this.control({
      'metabolitelist': {
	      select: this.onSelect,
	      deselect: this.onDeselect,
      },
      'metabolitelist button[action=clear]': {
        click: this.clearFilters
	    },
    });
  },
  onLoad: function(store) {
    console.log('Metabolite store loaded '+store.isLoaded);
    // TODO setChromatogramMarkersByMetaboliteFilter();
  },
  onSelect: function(rm, metabolite) {
    var me = this;
    var metid = metabolite.data.metid;
    console.log('Select metabolite '+metid);

    this.application.getController('Fragments').clearFragments();

// TODO
//     lc_chart.setLoading(true);
//     Ext.Ajax.request({
//       url: Ext.String.format(this.application.getUrls().extractedionchromatogram, metid),
//       success: function(response) {
//         lc_chart.setLoading(false);
//         var obj = Ext.decode(response.responseText);
//         metabolite.data.scans = obj.scans;
//         lc_chart.setExtractedIonChromatogram(obj.chromatogram);
//         if (obj.scans.length) {
//           // if one scan has already been selected test if its a member of the scans of selected metabolite
//           // if so then show fragments
//           if (
//             lc_chart.selectedscan != -1 &&
//             obj.scans.some(function(e) {
//               return (e.id == lc_chart.selectedscan);
//             })
//           ) {
//             console.log('Selected metabolite and its one scan is selected');
//             me.getController('Fragments').loadFragments(lc_chart.selectedscan, metid);
//           } else {
//             console.log('Selecting scans of metabolite');
//             if (lc_chart.selectedscan != -1) {
//               lc_chart.setMarkers(obj.scans);
//             } else {
//               var selectedScan = lc_chart.selectedscan;
//               lc_chart.setMarkers(obj.scans);
//               lc_chart.selectScans([selectedScan]);
//             }
//             // if metabolite has only one scan hit then show that scan and fragments
//             if (obj.scans.length == 1) {
//               var selectedScan = obj.scans[0].id;
//               // show scan where metabolite had its hit
//               console.log('show scan where metabolite had its hit');
//               lc_chart.selectScans([selectedScan]);
//               loadMSpectra1(selectedScan);
//               // show fragments with this metabolite and scan
//               console.log('show fragments of metabolite');
//               me.getController('Fragments').loadFragments(selectedScan, metid);
//             } else {
//               // multiple scans so mspectra1 should be empty
//               clearMSpectra1();
//             }
//           }
//         } else {
//           Ext.Msg.alert('Metabolite has no hits', 'The selected query/metabolite was not found in the ms data');
//           me.getMetaboliteList().getSelectionModel().deselectAll();
//         }
//       }
//     });
  },
  onDeselect: function(rm, r) {
    console.log('Deselect metabolite');
    // TODO
		//     setChromatogramMarkersByMetaboliteFilter();
		//     lc_chart.setExtractedIonChromatogram([]);
		this.application.getController('Fragments').clearFragments();
  },
  clearFilters: function() {
    console.log('Clear metabolite filter');
    this.getMetaboliteList().clearFilter();
    this.getMetabolitesStore().filter();
    // TODO
    //     lc_chart.setExtractedIonChromatogram([]);
    this.application.getController('Fragments').clearFragments();
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
});

Ext.define('Esc.msygma.view.fragment.Tree', {
  extend: 'Ext.tree.Panel',
  alias: 'widget.fragmenttree',
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
      id: 'fmolcol',
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
    this.getView().getPlugin('fmolcol').initCanvases();
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

    this.getFragmentsStore().on('load', this.onLoad);

    this.control({
      'fragmenttree': {
        select: this.onSelect,
        itemcollapse: this.onFragmentCollapse,
        itemexpand: this.onFragmentExpand
      }
    });
  },
  loadFragments: function (scanid, metid) {
    this.clearFragments();
    console.log('Show fragments of scan '+scanid+' metabolite '+metid);
    var store = this.getFragmentsStore();
    store.setProxy(this.fragmentProxyFactory(scanid, metid));
    store.load();
  },
  //need to change url of fragment proxy so use a factory to create a new proxy foreach scan/metabolite combo
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
    // TODO
		//     mspectras[1].setMarkers([]);
		//     for (var i = 2; i <= this.application.maxmslevel; i++) {
		//       clearMSpectra(i);
		//     }
  },
  onFragmentCollapse: function(fragment) {
    console.log('Collapsing fragment '+fragment.id);
    // on collapse clear child mspectra
// TODO
//     mspectras.forEach(function(ms, i) {
//       if (i > fragment.data.mslevel ) {
//         clearMSpectra(i);
//       }
//     });
  },
  onFragmentExpand: function(fragment) {
    console.log('Expanding fragment '+fragment.id);
    // on expand load child mspectra if needed
// TODO
    //     if (fragment.firstChild.data.scanid != mspectras[fragment.firstChild.data.mslevel].scanid) {
//       loadMSpectra(
//         fragment.firstChild.data.mslevel,
//         fragment.firstChild.data.scanid,
//         fragment.childNodes.map(function(r) { return {mz: r.data.mz}; })
//       );
//     }
  },
  onSelect: function(rm, r) {
    console.log('Selected fragment '+r.id);
    // select peak belonging to r
    // TODO
    //    mspectras[r.data.mslevel].selectPeak(r.data.mz);
    // show child mspectra of selected node or mz
//     if (!r.isLeaf()) {
//       // onselect then expand
//       if (!r.isExpanded()) {
//         r.expand();
//       } else {
//         this.onFragmentExpand(r);
//       }
//     }
//     // select peaks of parents of fragment in parent scans
//     if (r.data.mslevel==1) {
//       mspectras[1].selectPeak(r.data.mz);
//       for (var i = 2; i <= config.maxmslevel; i++) {
//         mspectras[i].clearPeakSelection();
//       }
//     } else if (r.data.mslevel==2) {
//       for (var i = 3; i <= config.maxmslevel; i++) {
//         mspectras[i].clearPeakSelection();
//       }
//       mspectras[2].selectPeak(r.data.mz);
//       mspectras[1].selectPeak(r.parentNode.data.mz);
//     } else if (r.data.mslevel>=3) {
//       // TODO make selecting parent Node work for mslevel>3
//       mspectras[2].selectPeak(r.parentNode.data.mz);
//       mspectras[1].selectPeak(r.parentNode.parentNode.data.mz);
//     }
  },
  onLoad: function(t, n, rs) {
    var me = this;
    // show peaks in lvl1 scan
    // TODO
//     if ('id' in n.data && n.data.id == 'root') {
//       console.log('Loaded metabolite as fragment');
//       // add mz of metabolites as markers to lvl1 scan
//       mspectras[1].setMarkers(
//         rs.map(function(r) { return {mz: r.data.mz}; })
//       );
//       var metabolite_fragment = rs[0];
//       mspectras[1].selectPeak(metabolite_fragment.data.mz);
//       if (metabolite_fragment.hasChildNodes()) {
//         loadMSpectra2(
//           metabolite_fragment.childNodes[0].data.scanid,
//           metabolite_fragment.childNodes.map(function(r) { return {mz: r.data.mz}; }),
//           function() {
//             // fgrid.refresh event is called before canvas have been rendered
//             // force fragment molecule rendering, hopyfully canvas have been rendered after spectra has been loaded
//             me.getFragmentTree().initMolecules();
//           }
//         );
//       }
//     } else if (n.data.mslevel == 1) {
//       console.log('Loaded lvl2 fragments of metabolite ');
//       // load the scan of first child
//       // add mz of metabolites as markers to lvl2 scan
//       loadMSpectra2(
//         rs[0].data.scanid,
//         rs.map(function(r) { return {mz: r.data.mz}; }),
//         function() {
//           // fgrid.refresh event is called before canvas have been rendered
//           // force fragment molecule rendering, hopyfully canvas have been rendered after spectra has been loaded
//           me.getFragmentTree().initMolecules();
//         }
//       );
//       mspectras[1].selectPeak(n.data.mz);
//     } else if (n.data.mslevel >= 2) {
//       console.log('Loaded lvl'+(n.data.mslevel+1)+' fragments of metabolite ');
//       // load the scan of first child
//       // add mz of metabolites as markers to lvl3 scan
//       loadMSpectra(
//         n.data.mslevel+1,
//         rs[0].data.scanid,
//         rs.map(function(r) { return {mz: r.data.mz}; }),
//         function() {
//           // fgrid.refresh event is called before canvas have been rendered
//           // force fragment molecule rendering, hopyfully canvas have been rendered after spectra has been loaded
//           me.getFragmentTree().initMolecules();
//         }
//       );
//       // TODO select parent peaks if n.data.mslevel>2
//       mspectras[2].selectPeak(n.data.mz);
//     }
  }
});

Ext.define('Esc.msygma.controller.Scans', {
  extend: 'Ext.app.Controller',
  init: function() {
    console.log('Scans controller init');
  }
});

Ext.define('Esc.msygma.resultsApp', { extend:'Ext.app.Application',
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
  applyMaxmslevel: function(val) {
    console.log('Apply maxmslevel');
    return val;
  },
  launch: function() {
    console.log('Launch app');
    var me = this;
    var config = me.config;

    function clearMSpectra(mslevel) {
      mspectras[mslevel].setData([]);
      mspectras[mslevel].scanid = -1;
      Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Scan ... (Level '+mslevel+')');
    }

    function loadChildMSpectraOfFragment(node) {
      // on expand load child mspectra if needed
      if (node.firstChild.data.scanid != mspectras[node.firstChild.data.mslevel].scanid) {
        loadMSpectra(
          node.firstChild.data.mslevel,
          node.firstChild.data.scanid,
          node.childNodes.map(function(r) { return {mz: r.data.mz}; })
        );
      }
    }

    /**
     * When user selects fragment in tree then select the peak in the mspectra
     */
    function selectPeakInMSpectra(rm, r) {
      console.log('Selected '+r.id);
      // select peak belonging to r
      mspectras[r.data.mslevel].selectPeak(r.data.mz);
      // show child mspectra of selected node or mz
      if (!r.isLeaf()) {
        // onselect then expand
        if (!r.isExpanded()) {
          r.expand();
        } else {
          loadChildMSpectraOfFragment(r);
        }
      }
      // select peaks of parents of fragment in parent scans
      if (r.data.mslevel==1) {
        mspectras[1].selectPeak(r.data.mz);
        for (var i = 2; i <= config.maxmslevel; i++) {
          mspectras[i].clearPeakSelection();
        }
      } else if (r.data.mslevel==2) {
        for (var i = 3; i <= config.maxmslevel; i++) {
          mspectras[i].clearPeakSelection();
        }
        mspectras[2].selectPeak(r.data.mz);
        mspectras[1].selectPeak(r.parentNode.data.mz);
      } else if (r.data.mslevel>=3) {
        // TODO make selecting parent Node work for mslevel>3
        mspectras[2].selectPeak(r.parentNode.data.mz);
        mspectras[1].selectPeak(r.parentNode.parentNode.data.mz);
      }
    }

    /**
     * When user selects peak in spectra then select the fragment belonging to peak in fragment tree
     */
    function selectFragmentInTree(mz, mslevel) {
      console.log('Selected peak in lvl'+mslevel+' mspectra with m/z = '+mz);
      // find fragment based on mz + mslevel
      var node = me.getController('Fragments').getFragmentsStore().getRootNode().findChildBy(function(n) {
        return (n.data.mslevel == mslevel && n.data.mz == mz);
      }, false, true);
      me.getController('Fragments').getFragmentTree().getSelectionModel().select([node]);
      if (!node.isLeaf()) {
        if (!node.isExpanded()) {
          node.expand();
        } else {
          loadChildMSpectraOfFragment(node);
        }
      } else {
        // clear product scans
      }
      // TODO unselect peaks of child scans
    }

    /**
     * Removes scan filter from metabolite store
     */
    function removeScanFilter() {
      var scanfilter;
      // see if already filtered on scanid then remove old filter
      if ('scanid' in me.getController('Metabolites').getMetabolitesStore().getProxy().extraParams) {
        delete(me.getController('Metabolites').getMetabolitesStore().getProxy().extraParams.scanid);
        me.getController('Metabolites').getMetabolitesStore().loadPage(1);
      }
    }

    function loadMSpectra1(scanid, onload) {
      loadMSpectra(1, scanid, [], onload);
    }

    function loadMSpectra2(scanid, markers, onload) {
      loadMSpectra(2, scanid, markers, onload);
    }

    function loadMSpectra(mslevel, scanid, markers, onload) {
      console.log('Loading msspectra level '+mslevel+' with id '+scanid);
      mspectras[mslevel].setLoading(true);
      d3.json(Ext.String.format(config.urls.mspectra, scanid, mslevel), function(data) {
        if (!data) {
          Ext.MessageBox.show({
            title: 'Unable to find scan',
            msg: 'Level '+mslevel+' scan with id '+scanid+' was not found',
            buttons: Ext.MessageBox.OK,
            icon: Ext.MessageBox.ERROR
          });
          mspectras[mslevel].setLoading(false);
          return;
        }
        Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Scan '+scanid+' (Level '+mslevel+')');
        mspectras[mslevel].setLoading(false);
        mspectras[mslevel].scanid = scanid;
        mspectras[mslevel].cutoff = data.cutoff;
        mspectras[mslevel].setData(data.peaks);
        mspectras[mslevel].setMarkers(markers);
        if (onload) {
          onload();
        }
      });
    }

    function selectScan(scanid) {
      console.log('select scan '+scanid);
      removeScanFilter();
      me.getController('Fragments').clearFragments();
      // if metabolite has been selected
      // and scanid is hit of metabolite then show fragments with this metabolite and scan
      // else filter mstore and load scan
      if (
          me.getController('Metabolites').getMetaboliteList().getSelectionModel().selected.getCount() > 0
          &&
          me.getController('Metabolites').getMetaboliteList().getSelectionModel().selected.getAt(0).data.scans.some(
            function(e) { return (e.id == scanid); }
          )
      ) {
        loadMSpectra1(scanid, function() {
          me.getController('Fragments').loadFragments(
              scanid,
              me.getController('Metabolites').getMetaboliteList().getSelectionModel().selected.getAt(0).data.metid
          );
        });
      } else {
        loadMSpectra1(scanid, function() {
          me.getController('Metabolites').getMetabolitesStore().getProxy().extraParams.scanid = scanid;
          me.getController('Metabolites').getMetabolitesStore().loadPage(1);
        });
        // TODO in spectra add markers for metabolites present in scan
      }
    }

    function clearMSpectra1() {
      clearMSpectra(1);
    }

    function unSelectScan() {
      removeScanFilter();
      me.getController('Fragments').clearFragments();
      clearMSpectra1();
    }

    function setChromatogramMarkersByMetaboliteFilter() {
      var store = me.getController('Metabolites').getMetabolitesStore();
      if (store.isLoaded && lc_chart.hasData()) {
        console.log('Setting chromatogram markers');
        var markers = store.getProxy().getReader().rawData.scans;
        lc_chart.setMarkers(markers);
      }
    }

    lc_chart = Ext.create('Ext.esc.Chromatogram', {
      cutoff: config.ms_intensity_cutoff,
      emptyText: 'Loading chromatogram ...',
      listeners: {
        selectscan: selectScan,
        unselectscan: unSelectScan
      }
    });
    lc_chart.setLoading(true);
    d3.json(config.urls.chromatogram, function(data) {
      lc_chart.setLoading(false);
      lc_chart.setData(data);
      setChromatogramMarkersByMetaboliteFilter();
    });

    var msspectrapanels = [];
    mspectras = [];
    for (var i = 1; i <= config.maxmslevel; i++) {
      mspectras[i] = Ext.create('Ext.esc.MSpectra', {
        emptyText: (
            i==1 ?
            'Select a scan in the chromatogram' :
            'Select a fragment to show its level '+i+' scan'
        ),
        listeners: {
          selectpeak: function(mz) {
            selectFragmentInTree(mz, i);
          }
        }
      });
      msspectrapanels.push({
        title: 'Scan ... (Level '+i+')',
        id: 'mspectra'+i+'panel',
        collapsible: true,
        items: mspectras[i]
      });
    }

    var master_side = Ext.create('Ext.panel.Panel', {
      // master side
      region: 'center',
      layout: 'border',
      border: false,
      items:[{
        region:'center',
        title: 'Query molecules & Metabolites',
        border: false,
        xtype: 'metabolitelist'
      },{
        title:'Chromatogram',
        region:'south',
        hideCollapseTool: true,
        collapsible: true,
        height: '50%',
        split: true,
        layout: 'fit',
        items:[lc_chart],
        border: false,
        tools:[{
          type:'search',
          tooltip: 'Select lvl1 scan by identifier',
          handler: function() {
            Ext.MessageBox.prompt('Scan#', 'Please enter a level 1 scan identifier:', function(b,v) {
              if (b != 'cancel' && v) {
               v = v*1;
               lc_chart.selectScans([v]);
               selectScan(v);
              }
            });
          }
        },{
          type:'refresh',
          tooltip:'Clear scan selection',
          handler: function() {
            lc_chart.clearScanSelection();
            unSelectScan();
          }
        }]
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
        title: 'Fragments',
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

app = Ext.create('Esc.msygma.resultsApp',{
	maxmslevel: ${maxmslevel},
	ms_intensity_cutoff: ${run.ms_intensity_cutoff},
	urls: {
		fragments: '${request.application_url}/fragments/{0}/{1}.json',
		mspectra: '${request.application_url}/mspectra/{0}.json?mslevel={1}',
		extractedionchromatogram: '${request.application_url}/extractedionchromatogram/{0}.json',
		metabolites: '${request.route_url('metabolites.json')}',
		chromatogram: '${request.route_url('chromatogram.json')}'
	}
});

</script>
</head>
<body>
</body>
</html>

