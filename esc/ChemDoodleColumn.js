// for post render clue see http://skirtlesden.com/ux/component-column
Ext.define('Ext.esc.ChemDoodleColumn', {
  extend: 'Ext.grid.column.Column',
  alias: ['widget.chemdoodlecolumn'],
  config: {
    canvasWidth:150,
    canvasHeight:100,
    // change when more then one chemdoodle columns in grid
    // used to find canvases belonging to this column
    canvasClass: 'x-chemdoodle-cols'
  },
  constructor: function(config) {
    this.callParent(arguments);
    // use chemdoodlecolumn subset of config to initconfig, otherwise problem with width not being able to be set
    this.initConfig({
      canvasWidth: config.canvasWidth || 150,
      canvasHeight: config.canvasHeight || 100,
      // change when more then one chemdoodle columns in grid
      // used to find canvases belonging to this column
      canvasClass: config.canvasClass || 'x-chemdoodle-cols'
    });
    return this;
  },
  renderer: function(v,meta,r,row,col,store,gridview) {
    // renderer returns html string
    // we want to run js on created dom, but that must be delayed until dom is complete.
    // create canvas with right dimensions using id unique to this gridview and cell (row/col)
    var c = this.columns[col];
    return '<canvas width='+c.getCanvasWidth()+' height='+c.getCanvasHeight()+' class="'+c.getCanvasClass()+'" id="'+gridview.id+'-'+row+'-'+col+'"></canvas>';
  },
  // use grid plugin so canvas can be selected with grid as root
  init: function(grid) {
    var self = this;
    this.grid = grid;
    // after canvas tag has been added to dom
    // find canvas tags and paint with chemdoodle
    this.grid.getView().on('refresh', function() {
      Ext.Array.forEach(
          Ext.DomQuery.select('canvas[class*="'+this.getCanvasClass()+'"]',this.grid.getView().dom),
          function(canvas) {
            var canvasid = canvas.id;
            var rowid = canvasid.replace(this.grid.getView().id+'-','');
            rowid = rowid.replace(/-\d$/,'');
            var row = this.grid.getStore().getAt(rowid);
            this.initCanvas(canvasid,this.getCanvasWidth(),this.getCanvasHeight(),row.data[this.dataIndex],row);
          }, 
          self
      );
    }, this);
  },
  initCanvas: function(id, width, height, value, record) {
    var c = new ChemDoodle.ViewerCanvas(id, width, height);
    c.loadMolecule(ChemDoodle.readMOL(value));
  } 
});