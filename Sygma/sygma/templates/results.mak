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
}

line.peak {
  stroke-width: 2px;
}

path.marker {
  stroke: lightgreen;
  fill: none;
  stroke-width: 1.5px;
}

.selectedscan, .selectedpeak {
  stroke: darkgreen !important;
  fill: darkgreen !important;
}

#fragmentgrid .x-grid-cell-inner{
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
Ext.require([
    'Ext.ux.grid.FiltersFeature',
]);

Ext.onReady(function () {

  function clearFragments() {
    console.log('Clearing fragments and mspectra');
    fstore.getRootNode().removeAll();
    mspectras[1].setMarkers([]);
    % for i in range(2,maxmslevel+1):
    clearMSpectra(${i});
    % endfor
  }

  function clearMSpectra(mslevel) {
    mspectras[mslevel].setData([]);
    mspectras[mslevel].scanid = -1;
    Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Scan ... (Level '+mslevel+')');
  }

  function loadFragments(scanid, metid) {
    clearFragments();
    console.log('Show fragments of scan '+scanid+' metabolite '+metid);
    fstore.setProxy(fragmentProxy(scanid, metid));
    fstore.load();
  }

  // need to change url of fragment proxy so use a factory to create a new proxy foreach scan/metabolite combo
  function fragmentProxy(scanid, metid) {
    return Ext.create('Ext.data.proxy.Ajax', {
      // url is build when scan and metabolite are selected
      url: '${request.application_url}/fragments/'+scanid+'/'+metid+'.json',
      reader: {
          type: 'json',
          root: 'children',
          idProperty: 'fragid'
      }
    });
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
      % for i in range(2,maxmslevel+1):
      mspectras[${i}].clearPeakSelection();
      % endfor
    } else if (r.data.mslevel==2) {
      % for i in range(3,maxmslevel+1):
        mspectras[${i}].clearPeakSelection();
      % endfor
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
    var node = fstore.getRootNode().findChildBy(function(n) {
      return (n.data.mslevel == mslevel && n.data.mz == mz);
    }, false, true);
    fgrid.getSelectionModel().select([node]);
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
    if ('scanid' in mstore.getProxy().extraParams) {
      delete(mstore.getProxy().extraParams.scanid);
      mstore.loadPage(1);
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
    d3.json('${request.application_url}/mspectra/'+scanid+'.json?mslevel='+mslevel, function(data) {
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
    clearFragments();
    // if metabolite has been selected
    // and scanid is hit of metabolite then show fragments with this metabolite and scan
    // else filter mstore and load scan
    if (
        mgrid.getSelectionModel().selected.getCount() > 0
        &&
        mgrid.getSelectionModel().selected.getAt(0).data.scans.some(
          function(e) { return (e.id == scanid); }
        )
    ) {
      loadMSpectra1(scanid, function() {
        loadFragments(scanid, mgrid.getSelectionModel().selected.getAt(0).data.metid);
      });
    } else {
      loadMSpectra1(scanid, function() {
        mstore.getProxy().extraParams.scanid = scanid;
        mstore.loadPage(1);
      });
      // TODO in spectra add markers for metabolites present in scan
    }
  }

  function clearMSpectra1() {
    clearMSpectra(1);
  }

  function unSelectScan() {
    removeScanFilter();
    clearFragments();
    clearMSpectra1();
  }

  function selectMetabolite(metabolite) {
    var metid = metabolite.data.metid;
    console.log('Selected metabolite '+metid);
    clearFragments();
    lc_chart.setLoading(true);
    Ext.Ajax.request({
      url: '${request.application_url}/extractedionchromatogram/'+metid+'.json',
      success: function(response) {
        lc_chart.setLoading(false);
        var obj = Ext.decode(response.responseText);
        metabolite.data.scans = obj.scans;
        lc_chart.setExtractedIonChromatogram(obj.chromatogram);
        if (obj.scans.length) {
          // if one scan has already been selected test if its a member of the scans of selected metabolite
          // if so then show fragments
          if (
            lc_chart.selectedscan != -1 &&
            obj.scans.some(function(e) {
              return (e.id == lc_chart.selectedscan);
            })
          ) {
            console.log('Selected metabolite and its one scan is selected');
            loadFragments(lc_chart.selectedscan, metid);
          } else {
            console.log('Selecting scans of metabolite');
            if (lc_chart.selectedscan != -1) {
              lc_chart.setMarkers(obj.scans);
            } else {
              var selectedScan = lc_chart.selectedscan;
              lc_chart.setMarkers(obj.scans);
              lc_chart.selectScans([selectedScan]);
            }
            // if metabolite has only one scan hit then show that scan and fragments
            if (obj.scans.length == 1) {
              var selectedScan = obj.scans[0].id;
              // show scan where metabolite had its hit
              console.log('show scan where metabolite had its hit');
              lc_chart.selectScans([selectedScan]);
              loadMSpectra1(selectedScan);
              // show fragments with this metabolite and scan
              console.log('show fragments of metabolite');
              loadFragments(selectedScan, metid);
            } else {
              // multiple scans so mspectra1 should be empty
              clearMSpectra1();
            }
          }
        } else {
          Ext.Msg.alert('Metabolite has no hits', 'The selected query/metabolite was not found in the ms data');
          mgrid.getSelectionModel().deselectAll();
        }
      }
    });
  }

  function setChromatogramMarkersByMetaboliteFilter() {
    var store = Ext.StoreMgr.get('metabolites');
    if (store.isLoaded && lc_chart.hasData()) {
      console.log('Setting chromatogram markers');
      var markers = store.getProxy().getReader().rawData.scans;
      lc_chart.setMarkers(markers);
    }
  }

  Ext.define('Metabolite', {
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

  Ext.define('Fragment', {
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

  var pageSize = 10;
  var mstore = Ext.create('Ext.data.Store', {
    storeId:'metabolites',
    model:'Metabolite',
    pageSize: pageSize,
    proxy: {
        type: 'ajax',
        url: '${request.route_url('metabolites.json')}',
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
    autoLoad: true,
    remoteSort: true,
    remoteFilter: true,
    isLoaded: false,
    listeners: {
        load: function(store) {
          this.isLoaded = true;
          setChromatogramMarkersByMetaboliteFilter();
        }
    }
 });

  var molcol = Ext.create('Ext.esc.ChemDoodleColumn', {
   text: 'Molecule', dataIndex: 'mol',
   width: 162
  });

  var mfilters = Ext.create('Ext.ux.grid.FiltersFeature',{
    encode: true,
  });

  var mselmodel = Ext.create('Ext.selection.CheckboxModel', {
    allowDeselect: true,
    mode: 'SINGLE'
  });

  var mgrid = Ext.create('Ext.grid.Panel', {
    id: 'metabolitegrid',
    store: mstore,
    selModel: mselmodel,
    columns: [
      {text: 'ID', dataIndex: 'metid', hidden: true},
      molcol,
      {text: 'Level', dataIndex: 'level', filter: { type: 'list',  options: [0,1,2,3] }, hidden:true},
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
    scroll: false,
    viewConfig: {
      autoScroll: true,
    },
    pageSize: pageSize,
    dockedItems: [{
      xtype: 'pagingtoolbar',
      store: mstore,   // same store GridPanel is using
      dock: 'bottom',
      displayInfo: true,
      items: [{
        text: 'Clear filters',
        handler: function() {
          mfilters.clearFilters();
          mstore.filter();
          lc_chart.setExtractedIonChromatogram([]);
          clearFragments();
        }
      }]
    }],
    plugins: [molcol],
    features: [mfilters],
    listeners: {
      deselect: function(m,rs) {
        if (rs.data.scans.length > 1) {
        }
        setChromatogramMarkersByMetaboliteFilter();
        lc_chart.setExtractedIonChromatogram([]);
        clearFragments();
      },
      select: function(rm,r) {
        selectMetabolite(r);
      }
    }
  });

  // atoms property is array filled with fragment atoms that need to be black
  // bonds having both atoms in array are black
  var hlspecs = new ChemDoodle.structures.VisualSpecifications();
  hlspecs.bonds_color = 'black';
  hlspecs.atoms_color = 'black';
  fmolcol = Ext.create('Ext.esc.ChemDoodleColumn', {
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

  var fstore = Ext.create('Ext.data.TreeStore', {
    storeId:'fragmentStore',
    model:'Fragment',
    proxy: fragmentProxy(),
    autoLoad: false,
    root: { children : [] }, // prevent tree from autoloading
    listeners: {
      load: function(t, n, rs) {
        // show peaks in lvl1 scan
        if ('id' in n.data && n.data.id == 'root') {
          console.log('Loaded metabolite as fragment');
          // add mz of metabolites as markers to lvl1 scan
          mspectras[1].setMarkers(
            rs.map(function(r) { return {mz: r.data.mz}; })
          );
          var metabolite_fragment = rs[0];
          mspectras[1].selectPeak(metabolite_fragment.data.mz);
          if (metabolite_fragment.hasChildNodes()) {
            loadMSpectra2(
              metabolite_fragment.childNodes[0].data.scanid,
              metabolite_fragment.childNodes.map(function(r) { return {mz: r.data.mz}; }),
              function() {
                // fgrid.refresh event is called before canvas have been rendered
                // force fragment molecule rendering, hopyfully canvas have been rendered after spectra has been loaded
                fmolcol.initCanvases();
              }
            );
          }
        } else if (n.data.mslevel == 1) {
          console.log('Loaded lvl2 fragments of metabolite ');
          // load the scan of first child
          // add mz of metabolites as markers to lvl2 scan
          loadMSpectra2(
            rs[0].data.scanid,
            rs.map(function(r) { return {mz: r.data.mz}; }),
            function() {
              // fgrid.refresh event is called before canvas have been rendered
              // force fragment molecule rendering, hopyfully canvas have been rendered after spectra has been loaded
              fmolcol.initCanvases();
            }
          );
          mspectras[1].selectPeak(n.data.mz);
        } else if (n.data.mslevel >= 2) {
          console.log('Loaded lvl'+(n.data.mslevel+1)+' fragments of metabolite ');
          // load the scan of first child
          // add mz of metabolites as markers to lvl3 scan
          loadMSpectra(
            n.data.mslevel+1,
            rs[0].data.scanid,
            rs.map(function(r) { return {mz: r.data.mz}; }),
            function() {
              // fgrid.refresh event is called before canvas have been rendered
              // force fragment molecule rendering, hopyfully canvas have been rendered after spectra has been loaded
              fmolcol.initCanvases();
            }
          );
          // TODO select parent peaks if n.data.mslevel>2
          mspectras[2].selectPeak(n.data.mz);
        }
      }
    },
    // TreeStore and Store have different function to fetch record by id, add getById to TreeStore
    getById: function(id) {
      return this.getNodeById(id);
    }
  });

  var fgrid = Ext.create('Ext.tree.Panel', {
    id: 'fragmentgrid',
    store: fstore,
    selType:'checkboxmodel',
    multiSelect: false,
    rootVisible: false,
    singleExpand: true,
    scroll: false,
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
    viewConfig: {
      // animate is default true causing refresh event to be blocked
      // we use refresh event to render molecules
      // so after expanding a node the refresh was not fired causing all prev. rendered mols to disappear
      // now we turn off animate, so refresh events are fired and mols can be rendered
      // But still some nodes have no molecule and be fixes by collapse/expand or fmolcol.initCanvases()
      animate: false,
      autoScroll: true,
      blockRefresh: false,
      emptyText: 'Select a metabolite and scan, to show its fragments',
    },
    listeners: {
      select: selectPeakInMSpectra,
      itemcollapse: function(node) {
        // on collapse clear child mspectra
        mspectras.forEach(function(ms, i) {
          if (i > node.data.mslevel ) {
            clearMSpectra(i);
          }
        });
      },
      itemexpand: loadChildMSpectraOfFragment
    }
  });

  lc_chart = Ext.create('Ext.esc.Chromatogram', {
    cutoff: ${run.ms_intensity_cutoff},
    emptyText: 'Loading chromatogram ...',
    listeners: {
      selectscan: selectScan,
      unselectscan: unSelectScan
    }
  });
  lc_chart.setLoading(true);
  d3.json('${request.route_url('chromatogram.json')}', function(data) {
    lc_chart.setLoading(false);
    lc_chart.setData(data);
    setChromatogramMarkersByMetaboliteFilter();
  });

  var msspectrapanels = [];
  mspectras = [];
  % for i in range(1,maxmslevel+1):
  mspectras[${i}] = Ext.create('Ext.esc.MSpectra', {
    emptyText:
      % if i==1:
      'Select a scan in the chromatogram',
      % else:
      'Select a fragment to show its level ${i} scan',
      % endif
    listeners: {
      selectpeak: function(mz) {
        selectFragmentInTree(mz, ${i});
      }
    }
  });
  msspectrapanels.push({
    title: 'Scan ... (Level ${i})',
    id: 'mspectra${i}panel',
    collapsible: true,
    items: mspectras[${i}]
  });
  % endfor

  var master_side = Ext.create('Ext.panel.Panel', {
    // master side
    region: 'center',
    layout: 'border',
    border: false,
    items:[{
      region:'center',
      title: 'Query molecules & Metabolites',
      layout: 'fit',
      border: false,
      items: [
        mgrid,
      ]
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
      id: 'fragmentspanel',
      title: 'Fragments',
      items: [fgrid],
      layout: 'fit',
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

  Ext.QuickTips.init();
});
</script>
</head>
<body>
</body>
</html>

