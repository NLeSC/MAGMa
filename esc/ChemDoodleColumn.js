/**
 * @class Ext.esc.ChemDoodleColumn
 * @extends Ext.grid.column.Column
 * 
 * <p>A Column definition class which renders a chemdoodle canvas</p>
 * 
 * ## Code
 *     Ext.create('Ext.data.Store', {
 *        storeId:'sampleStore',
 *        fields:[
 *            {name: 'name', type: 'string'},
 *            {name: 'mol', type: 'string'}
 *        ],
 *        data:{'items':[
 *            name: 'Cafeine',
 *            mol: 'Molecule Name\n  CHEMDOOD08070920033D 0   0.00000     0.00000     0\n[Insert Comment Here]\n 14 15  0  0  0  0  0  0  0  0  1 V2000\n   -0.3318    2.0000    0.0000   O 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.3318    1.0000    0.0000   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -1.1980    0.5000    0.0000   N 0  0  0  1  0  0  0  0  0  0  0  0\n    0.5342    0.5000    0.0000   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -1.1980   -0.5000    0.0000   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -2.0640    1.0000    0.0000   C 0  0  0  4  0  0  0  0  0  0  0  0\n    1.4804    0.8047    0.0000   N 0  0  0  1  0  0  0  0  0  0  0  0\n    0.5342   -0.5000    0.0000   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -2.0640   -1.0000    0.0000   O 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.3318   -1.0000    0.0000   N 0  0  0  1  0  0  0  0  0  0  0  0\n    2.0640   -0.0000    0.0000   C 0  0  0  2  0  0  0  0  0  0  0  0\n    1.7910    1.7553    0.0000   C 0  0  0  4  0  0  0  0  0  0  0  0\n    1.4804   -0.8047    0.0000   N 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.3318   -2.0000    0.0000   C 0  0  0  4  0  0  0  0  0  0  0  0\n  1  2  2  0  0  0  0\n  3  2  1  0  0  0  0\n  4  2  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  6  1  0  0  0  0\n  7  4  1  0  0  0  0\n  4  8  2  0  0  0  0\n  9  5  2  0  0  0  0\n 10  5  1  0  0  0  0\n 10  8  1  0  0  0  0\n  7 11  1  0  0  0  0\n  7 12  1  0  0  0  0\n 13  8  1  0  0  0  0\n 13 11  2  0  0  0  0\n 10 14  1  0  0  0  0\nM  END\n> <DATE>\n07-08-2009\n'
 *        },{
 *            name: 'Morphine',
 *            mol: 'Molecule Name\n  CHEMDOOD01231009093D 0   0.00000     0.00000     0\n[Insert Comment Here]\n 40 44  0  0  0  0  0  0  0  0  1 V2000\n    1.7910    1.5052    1.7843   C 0  0  0  1  0  0  0  0  0  0  0  0\n    1.8814   0.4102    1.9755   C 0  0  0  1  0  0  0  0  0  0  0  0\n    1.2267    1.8815    0.8827   C 0  0  0  1  0  0  0  0  0  0  0  0\n    2.1840    2.0692    2.3166   H 0  0  0  1  0  0  0  0  0  0  0  0\n    2.4684    0.0562    2.8203   O 0  0  0  1  0  0  0  0  0  0  0 0\n    1.3897   -0.2898    1.2829   C 0  0  0  1  0  0  0  0  0  0  0  0\n    1.2058    2.7360    0.7207   H 0  0  0  1  0  0  0  0 0  0  0  0\n    0.7230    1.1614    0.1905   C 0  0  0  1  0  0  0  0  0  0  0  0\n    2.3966   -0.7208    2.8246   H 0  0  0  1  0 0  0  0  0  0  0  0\n    0.7850    0.0855    0.4453   C 0  0  0  1  0  0  0  0  0  0  0  0\n    1.4809   -1.3816    1.2883   O 0  0 0  1  0  0  0  0  0  0  0  0\n    0.2152    1.4793   -0.8555   C 0  0  0  1  0  0  0  0  0  0  0  0\n    0.3030   -0.8185   -0.1689   C 0  0  0  1  0  0  0  0  0  0  0  0\n    1.0589   -1.7109    0.2315   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.3585    2.1265   -0.7192   H 0  0  0  1  0  0  0  0  0  0  0  0\n    0.8523    1.8147   -1.3591   H 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.3072    0.5314   -1.4938   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.8642   -0.9923    0.1890   C 0  0  0  1  0  0  0  0  0  0  0  0\n    0.3287   -0.5398   -1.3656   C 0  0  0  1  0  0  0  0  0  0  0  0\n    0.6138   -2.4534    0.3781   H 0  0  0  1  0  0  0  0  0  0 0  0\n    1.9958   -1.9810   -0.5446   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -1.4413    0.2916   -1.2228   N 0  0  0  1  0  0  0 0  0  0  0  0\n   -0.2905    0.7503   -2.3452   H 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.9165   -1.1267    1.0554   H 0  0  0  1 0  0  0  0  0  0  0  0\n   -1.5945   -0.0502   -0.1099   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -1.1717   -1.7202   -0.1951   H 0 0  0  1  0  0  0  0  0  0  0  0\n    1.4542   -0.5270   -1.8056   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.0941   -1.1628   -1.8222   H 0  0  0  1  0  0  0  0  0  0  0  0\n    1.8036   -2.7360   -0.9489   H 0  0  0  1  0  0  0  0  0  0  0  0\n    2.2104   -1.1828   -1.4164   C 0  0  0  1  0  0  0  0  0  0  0  0\n    2.9720   -2.1604    0.0172   O 0  0  0  1  0  0  0  0  0  0  0  0\n   -2.1740    1.1577   -1.5022   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -2.4251   -0.3100    0.0050   H 0  0  0  1  0  0  0  0  0  0 0  0\n   -1.4674    0.6168    0.4458   H 0  0  0  1  0  0  0  0  0  0  0  0\n    1.6378   -0.0028   -2.4745   H 0  0  0  1  0  0  0 0  0  0  0  0\n    3.0043   -1.1859   -1.7720   H 0  0  0  1  0  0  0  0  0  0  0  0\n    2.9664   -1.6654    0.6201   H 0  0  0  1 0  0  0  0  0  0  0  0\n   -2.0760    1.8646   -0.9948   H 0  0  0  1  0  0  0  0  0  0  0  0\n   -3.0043    0.8869   -1.4262   H 0 0  0  1  0  0  0  0  0  0  0  0\n   -2.0673    1.3986   -2.3381   H 0  0  0  1  0  0  0  0  0  0  0  0\n  1  2  2  0  0  0  0\n  1  3  1 0  0  0  0\n  1  4  1  0  0  0  0\n  2  5  1  0  0  0  0\n  2  6  1  0  0  0  0\n  3  7  1  0  0  0  0\n  8  3  2  0  0  0  0\n  5  9  1  0 0  0  0\n  6 10  2  0  0  0  0\n  6 11  1  0  0  0  0\n 10  8  1  0  0  0  0\n  8 12  1  0  0  0  0\n 10 13  1  0  0  0  0\n 14 11  1 0  0  0  0\n 12 15  1  0  0  0  0\n 12 16  1  0  0  0  0\n 17 12  1  0  0  0  0\n 13 18  1  0  0  0  0\n 13 14  1  0  0  0  0\n 13 19  1  0  0  0  0\n 14 20  1  0  0  0  0\n 21 14  1  0  0  0  0\n 17 22  1  0  0  0  0\n 17 23  1  0  0  0  0\n 19 17  1  0  0  0 0\n 18 24  1  0  0  0  0\n 18 25  1  0  0  0  0\n 18 26  1  0  0  0  0\n 19 27  1  0  0  0  0\n 19 28  1  0  0  0  0\n 21 29  1  0 0  0  0\n 30 21  1  0  0  0  0\n 21 31  1  0  0  0  0\n 22 25  1  0  0  0  0\n 22 32  1  0  0  0  0\n 25 33  1  0  0  0  0\n 25 34 1  0  0  0  0\n 27 35  1  0  0  0  0\n 27 30  2  0  0  0  0\n 30 36  1  0  0  0  0\n 31 37  1  0  0  0  0\n 32 38  1  0  0  0  0\n 32 39  1  0  0  0  0\n 32 40  1  0  0  0  0\nM  END'
 *        ]},
 *        proxy: {
 *            type: 'memory',
 *            reader: {
 *                type: 'json',
 *                root: 'items'
 *            }
 *        }
 *     });
 *     
 *     var molcol = Ext.create('Ext.esc.ChemDoodleColumn', {
 *       text: 'Molecule', dataIndex: 'mol',
 *       width: 150
 *     });
 *     
 *     Ext.create('Ext.grid.Panel', {
 *         title: 'Chemdoodle Column Demo',
 *         store: Ext.data.StoreManager.lookup('sampleStore'),
 *         columns: [
 *             {text: 'Name',  dataIndex: 'name', flex: 1},
 *             molcol
 *         ],
 *         height: 200,
 *         width: 400,
 *         plugins:[molcol],
 *         renderTo: Ext.getBody()
 *     });
 *
 */
Ext.define('Ext.esc.ChemDoodleColumn', {
  extend: 'Ext.grid.column.Column',
  alias: ['widget.chemdoodlecolumn'],
  config: {
    /**
     * @cfg {Number} canvasWidth
     * The width of the chemdoodle canvas in pixels.
     */
    canvasWidth:150,
    /**
     * @cfg {Number} canvasHeight
     * The height of the chemdoodle canvas in pixels.
     */
    canvasHeight:100,
    /**
     * @cfg {String} canvasClass
     * Class of the chemdoodle canvas tag.
     * 
     *  Change when more then one different chemdoodle columns in grid
     */
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
  // for post render clue see http://skirtlesden.com/ux/component-column
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
  /**
   * @method
   * <p>Initializes a ChemDoodle.ViewerCanvas in the cell and loads 'value' as a molfile. 
   * Can be overwritten to use a different canvas, perform highlighting etc. Example</p>
   * 
   * <pre><code>{
    initCanvas: function(id, width, height, value, record){
        var c = new ChemDoodle.TransformCanvas(id, width, height,true);
        c.specs.bonds_color = 'red';
        c.loadMolecule(ChemDoodle.readMOL(value));
    }
}
   * </code></pre>
   * 
   * @param {String} id Canvas identifier  
   * @param {Number} width
   * @param {Number} height
   * @param {Mixed} value The data value for the current cell
   * @param {Ext.data.Model} record The record for the current row
   */
  initCanvas: function(id, width, height, value, record) {
    var c = new ChemDoodle.ViewerCanvas(id, width, height);
    c.loadMolecule(ChemDoodle.readMOL(value));
  } 
});