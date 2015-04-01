/**
 * Grid of molecules.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.molecule.List', {
  extend: 'Ext.grid.Panel',
  alias: 'widget.moleculelist',
  requires: [
    'Ext.ux.grid.FiltersFeature', 'Esc.chemdoodle.Column',
    'Ext.toolbar.Paging', 'Ext.grid.column.Boolean',
    'Ext.grid.column.Action', 'Ext.selection.CheckboxModel',
    'Ext.grid.column.Number',
    'Esc.magmaweb.view.molecule.ReactionColumn',
    'Esc.magmaweb.view.molecule.ReactionFilter'
  ],
  store: 'Molecules',
  viewConfig: {
    emptyText: 'No structures available: Add structures or relax filters'
  },
  selType: 'checkboxmodel',
  selModel: {
    allowDeselect: true,
    mode: 'SINGLE',
    showHeaderCheckbox: false
  },
  dockedItems: [{
    xtype: 'pagingtoolbar',
    store: 'Molecules',   // same store GridPanel is using
    dock: 'bottom',
    displayInfo: true,
    items: [{
      xtype: 'combo',
      width: 170,
      store: [
        [10, '10 molecules per page'],
        [25, '25 molecules per page'],
        [50, '50 molecules per page'],
        [100, '100 molecules per page'],
        [250, '250 molecules per page'],
        [500, '500 molecules per page'],
        [1000, '1000 molecules per page']
      ],
      forceSelection: true,
      triggerAction: 'all',
      action: 'pagesizeCombo'
    }]
  }],
  initComponent: function() {
    Ext.log({}, 'Init met grid');
    var me = this;
    var molcol = Ext.create('Esc.chemdoodle.Column', {
      text: 'Molecule', dataIndex: 'mol',
      width: 162,
      initCanvas: function(id, width, height, value, record) {
        var c = new ChemDoodle.ViewerCanvas(id, width, height);
        c.loadMolecule(ChemDoodle.readMOL(value));
        var tip = Ext.create('Ext.tip.ToolTip', {
          target: id,
          dismissDelay: 0,
          html: null,
          listeners: {
            show: function(tip) {
              tip.setSize(300, 300);
            },
            render: function(tip) {
              var c = new ChemDoodle.ViewerCanvas(id+'-'+tip.id, 290, 290);
              c.loadMolecule(ChemDoodle.readMOL(value));
             }
          }
        });
        tip.update('<canvas id="'+id+'-'+tip.id+'"></canvas>');
      }
    });

    var mfilters = Ext.create('Ext.ux.grid.FiltersFeature',{
      id: 'mfilter',
      encode: true
    });

    this.addEvents([
        /**
         * @event metabolize
         * Fired after metabolize column action is clicked
         * @param {Ext.data.Model} record The record to metabolize
         */
        'metabolize'
    ]);

    var numberFilter = {
        type: 'numeric', fieldCfg: {lt: {decimalPrecision: 18}, gt: {decimalPrecision: 18}, eq: {decimalPrecision: 18}}
    };
    Ext.apply(this, {
      columns: [
        {text: 'ID', width:40, dataIndex: 'molid', hidden: true, filter: { type: 'numeric' }},
        {text: 'Scans', width:50, dataIndex: 'nhits', filter: {
            type: 'numeric', value:{gt:0}, active: true
        }},
        {text: 'Assigned', width:60, dataIndex: 'assigned', hidden: false, xtype:'booleancolumn', trueText:'Yes', falseText:'No', filter: { type: 'boolean' }},
        {text: 'Candidate score', dataIndex: 'score', hidden: true, filter: numberFilter, xtype: 'numbercolumn', format: '0.00000'},
        molcol,
        {text: 'Inchikey', dataIndex: 'inchikey14', hidden:true},
        {text: 'Smiles', dataIndex: 'smiles', hidden:true, filter: { type: 'string' }},
        {text: 'Formula', width:100, dataIndex: 'formula', filter: { type: 'string' }},
        {text: 'Mass', width:80, dataIndex: 'mim', filter: { type: 'numeric' }, hidden: false, xtype: 'numbercolumn', format: '0.00000'},
        {text: '&Delta;Mass (ppm)', width:80, dataIndex: 'deltappm', hidden: true, filter: { type: 'numeric' }, xtype: 'numbercolumn', format: '0.00000'},
        {text: 'M/z', width:80, dataIndex: 'mz', hidden: true, filter: numberFilter , xtype: 'numbercolumn', format: '0.00000'},
        {text: 'Name', dataIndex: 'name', flex:1, filter: { type: 'string' }},
        {
            text: 'Reactions', dataIndex: 'reactionsequence', flex:1, filter: { type: 'reaction' },
            xtype: 'reactioncolumn'
        },
        {text: 'Refscore', width:80, dataIndex: 'refscore', filter: { type: 'numeric' }, xtype: 'numbercolumn', format: '0.00000'},
        {text: 'LogP', dataIndex: 'logp', filter: { type: 'numeric' }, hidden: true, xtype: 'numbercolumn', format: '0.00000'},
        {text: 'Reference', dataIndex: 'reference', filter: { type: 'string' }, sortable: false },
        {text: 'Predicted', dataIndex: 'predicted', xtype:'booleancolumn', hidden: true, trueText:'Yes', falseText:'No', filter: { type: 'boolean' }},
        {xtype: 'actioncolumn', width:30, text:'Commands', hidden: true,
            items: [{
                tooltip: 'Metabolize',
                iconCls: 'metabolize-col',
                handler: function(grid, rowIndex) {
                    me.fireEvent('metabolize', me.getStore().getAt(rowIndex));
                }
            }]
        }
      ],
      plugins: [molcol],
      features: [mfilters]
    });
    this.callParent(arguments);
  },
  /**
   * @return {Ext.ux.grid.FiltersFeature}
   * @private
   */
  getFilter: function() {
      return this.getView().getFeature('mfilter');
  },
  /**
   * Clears all filters applied to molecules
   */
  clearFilters: function() {
    this.getFilter().clearFilters();
  },
  /**
   * Get column with scores.
   *
   * Column should only be shown when a scan has been selected.
   *
   * @return {Ext.grid.column.Column}
   */
  getFragmentScoreColumn: function() {
      return this.columns.filter(function(c) { return (c.dataIndex == "score");})[0];
  },
  /**
   * Get column with deltappm.
   *
   * Column should only be shown when a scan has been selected.
   *
   * @return {Ext.grid.column.Column}
   */
  getFragmentDeltaPpmColumn: function() {
      return this.columns.filter(function(c) { return (c.dataIndex == "deltappm");})[0];
  },
  /**
   * Get column with commands.
   *
   * @return {Ext.grid.column.Action}
   */
  getCommandsColumn: function() {
      return this.columns.filter(function(c) { return (c.text == "Commands");})[0];
  },
  /**
   * Hides column with scores/deltappm and removes sort in header column.
   */
  hideFragmentScoreColumn: function() {
    var sortCls = [
        Ext.baseCSSPrefix + 'column-header-sort-ASC',
        Ext.baseCSSPrefix + 'column-header-sort-DESC'
    ];
    var score = this.getFragmentScoreColumn();
    score.hide();
    score.removeCls(sortCls);
    var deltappm = this.getFragmentDeltaPpmColumn()
    deltappm.hide();
    deltappm.removeCls(sortCls);
  },
  /**
   * Shows column with scores/deltappm.
   */
  showFragmentScoreColumn: function() {
    this.getFragmentScoreColumn().show();
    this.getFragmentDeltaPpmColumn().show();
  },
  /**
   * JSON encoded the filter data as query
   * @return {String}
   */
  getFilterQuery: function() {
      var filter = this.getFilter();
      return filter.buildQuery(filter.getFilterData());
  },
  /**
   * Array of dataindexes currently visible. In order they appear.
   * @return {Array}
   */
  getVisiblColumnIndices: function() {
      return this.getView().getHeaderCt().getVisibleGridColumns().map(function(v) {return v.dataIndex}).filter(function(v) { return v});
  },
  hideCommandsColumn: function() {
      var cmds = this.getCommandsColumn();
      cmds.disableAction(0);
      cmds.disable();
      cmds.hide();
  },
  setPageSize: function(size) {
      Ext.ComponentQuery.query('component[action=pagesizeCombo]')[0].select(size);
  },
  getMzFilter: function() {
    var filter = this.getFilter();
    return filter.filters.get('mz');
  },
  setMzFilterToEqual: function(mz) {
    var mzfilter = this.getMzFilter();
    if (!mzfilter.active) {
    	mzfilter.setActive(true);
    }
    mzfilter.setValue({'eq': mz});
  },
  clearMzFilter: function() {
    var mzfilter = this.getMzFilter();
    if (mzfilter.active) {
    	mzfilter.setActive(false);
    }
  }
});
