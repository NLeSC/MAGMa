/**
 * Fragment tree.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.fragment.Tree', {
  extend: 'Ext.tree.Panel',
  alias: 'widget.fragmenttree',
  title: 'Fragments',
  store: 'Fragments',
  selType: 'checkboxmodel',
  selModel: {
     allowDeselect: true,
     mode: 'SINGLE',
     showHeaderCheckbox: false
  },
  cls: 'fragmenttree', // So height of column can be set with .fragmenttree .x-grid-cell-inner{height: 106px !important;}
  multiSelect: false,
  rootVisible: false,
  singleExpand: true,
  rowLines: true,
  viewConfig: {
	getRowClass: function(record) {
		// Make transition between mslevel visible by giving even/odd mslevel different bg color
    	return record.get("mslevel") % 2 === 0 ? this.altRowCls : "";
	}
  },
  tools: [{
     type: 'save',
     disabled: true,
     tooltip: 'Save fragment tree'
  }],
  dockedItems: [{
      xtype: 'toolbar',
      dock: 'bottom',
      layout: {
        type: 'hbox',
        align: 'middle',
        pack: 'center'
      },
      items: [{
          iconCls: 'icon-connect',
          text: 'Assign molecule to peak',
          action: 'assign_struct2peak',
          disabled: true,
          enableToggle: true
      }]
  }],
  requires: [ 'Esc.chemdoodle.Column', 'Ext.selection.CheckboxModel' ],
  viewConfig: {
    // animate is default true causing refresh event to be blocked
    // we use refresh event to render molecules
    // so after expanding a node the refresh was not fired causing all prev. rendered mols to disappear
    // now we turn off animate, so refresh events are fired and mols can be rendered
    animate: false,
    blockRefresh: false,
    deferEmptyText: false,
    emptyText: 'Select a metabolite and scan, to show its fragments'
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
        var c = new ChemDoodle.ViewerCanvas(id, width, height);
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
        var tip = Ext.create('Ext.tip.ToolTip', {
            target: id,
            listeners: {
                render: function(tip,e) {
                    console.log('Drawing tooltip for '+id);
                    var c = new ChemDoodle.ViewerCanvas(id+'-'+tip.id, 300, 300);
                    c.specs.bonds_color = 'cyan';
                    c.specs.atoms_color = 'cyan';
                    c.loadMolecule(m);
                }
            }
         });
        // if canvas doesn't have tip.id then tooltip is rendered twice, once with molecule, once halfheight blue
        tip.update('<canvas id="'+id+'-'+tip.id+'"></canvas>');
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
        { text: 'MS Level', dataIndex: 'mslevel'},
        { text: 'Fragment atoms', dataIndex: 'atoms', hidden: true},
        { text: '&Delta;H', dataIndex: 'deltah'},
        { text: '&Delta;Mass (ppm)', dataIndex: 'deltappm', hidden: true}
      ],
      plugins: [fmolcol]
    });

    this.callParent(arguments);
  },
  /**
   * Forces all molecules on fragment tree to be drawn
   */
  initMolecules: function() {
    // force fragment molecule rendering, hopyfully canvas have been rendered after spectra has been loaded
    return this.getPlugin('fmolcol').initCanvases();
  },
  getAssignStruct2PeakButton: function() {
      return this.query('component[action=assign_struct2peak]')[0];
  }
});
