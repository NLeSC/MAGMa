describe('Esc.ChemDoodleColumn', function() {

    var caffeineMolFile = 'Molecule Name\n  CHEMDOOD08070920033D 0   0.00000     0.00000     0\n[Insert Comment Here]\n 14 15  0  0  0  0  0  0  0  0  1 V2000\n   -0.3318    2.0000    0.0000   O 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.3318    1.0000    0.0000   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -1.1980    0.5000    0.0000   N 0  0  0  1  0  0  0  0  0  0  0  0\n    0.5342    0.5000    0.0000   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -1.1980   -0.5000    0.0000   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -2.0640    1.0000    0.0000   C 0  0  0  4  0  0  0  0  0  0  0  0\n    1.4804    0.8047    0.0000   N 0  0  0  1  0  0  0  0  0  0  0  0\n    0.5342   -0.5000    0.0000   C 0  0  0  1  0  0  0  0  0  0  0  0\n   -2.0640   -1.0000    0.0000   O 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.3318   -1.0000    0.0000   N 0  0  0  1  0  0  0  0  0  0  0  0\n    2.0640   -0.0000    0.0000   C 0  0  0  2  0  0  0  0  0  0  0  0\n    1.7910    1.7553    0.0000   C 0  0  0  4  0  0  0  0  0  0  0  0\n    1.4804   -0.8047    0.0000   N 0  0  0  1  0  0  0  0  0  0  0  0\n   -0.3318   -2.0000    0.0000   C 0  0  0  4  0  0  0  0  0  0  0  0\n  1  2  2  0  0  0  0\n  3  2  1  0  0  0  0\n  4  2  1  0  0  0  0\n  3  5  1  0  0  0  0\n  3  6  1  0  0  0  0\n  7  4  1  0  0  0  0\n  4  8  2  0  0  0  0\n  9  5  2  0  0  0  0\n 10  5  1  0  0  0  0\n 10  8  1  0  0  0  0\n  7 11  1  0  0  0  0\n  7 12  1  0  0  0  0\n 13  8  1  0  0  0  0\n 13 11  2  0  0  0  0\n 10 14  1  0  0  0  0\nM  END\n> <DATE>\n07-08-2009\n';

    function createDefault() {
        return Ext.create('Esc.chemdoodle.Column',{
            text: 'Molecule', dataIndex: 'mol',
            width: 162
        });
    }

    describe('create', function() {
        it('default', function() {
            var col = createDefault();
            expect(col.canvasWidth).toEqual(150);
            expect(col.canvasHeight).toEqual(100);
            expect(col.canvasClass).toEqual('x-chemdoodle-cols');
        });

        it('custom', function() {
            var col = Ext.create('Esc.chemdoodle.Column',{
                text: 'Molecule', dataIndex: 'mol',
                width: 162, canvasWidth: 250, canvasHeight:200,
                canvasClass: 'x-chemdoodle-cols2'
            });
            expect(col.canvasWidth).toEqual(250);
            expect(col.canvasHeight).toEqual(200);
            expect(col.canvasClass).toEqual('x-chemdoodle-cols2');
        });
    });

    it('renderer', function() {
        var record = { mol:caffeineMolFile, internalId: 5 };
        var col = createDefault();
        // col.renderer is called by Ext.grid.Panel so mock a grid
        var grid = function() {};
        grid.columns = [];
        grid.columns[2] = col;
        var f = col.renderer;
        var rendered = f.call(col, record.mol, {}, record, 1, 2, {}, { id: 'gridview-1234'});
        expect(rendered).toEqual('<canvas width=150 height=100 class="x-chemdoodle-cols" id="gridview-1234-5-2"></canvas>');
    });

    it('init', function() {
        // inits so initCanvases is called when grid view fires refresh event
        var view = function() {};
        view.on = function() {};
        spyOn(view, 'on');
        var grid = function() {};
        grid.on = function() {};
        spyOn(grid, 'on');
        grid.getView = function() {
            return view;
        };
        var col = createDefault();
        col.init(grid);
        expect(view.on).toHaveBeenCalledWith('refresh',col.initCanvases, col);
        expect(grid.on).toHaveBeenCalledWith('afteritemcollapse',col.initCanvases, col);
        expect(col.grid).toEqual(grid);
    });

    describe('initCanvases', function() {
        it('nocanvases', function() {
            var col = createDefault();
            col.grid = function() { };
            col.grid.getView = function() {};
            // calls initCanvas on all canvas tags with class="x-chemdoodle-cols"
            var canvases = [];
            spyOn(col.grid,'getView').andReturn(function() { dom = 'dom' });
            spyOn(Ext.DomQuery,'select').andReturn(canvases);
            col.initCanvases();
            expect(Ext.DomQuery.select).toHaveBeenCalled();
        });

        it('1canvas', function() {
            // calls initCanvas on all canvas tags with class="x-chemdoodle-cols"
            // canvas id contains a record id
            var view = function() {};
            view.dom = function() { return 'dom'; };
            view.id = 'gridview-1234';
            var store = function() {};
            var record = { data:{ 'mol': caffeineMolFile } };
            store.data = { getByKey: function() { return record } };
            var col = createDefault();
            col.grid = function() { };
            col.grid.getView = function() { return view };
            col.grid.getStore = function() { return store };
            var canvases = [{ id:'gridview-1234-5-2' }];

            spyOn(Ext.DomQuery,'select').andReturn(canvases);
            spyOn(col,'initCanvas');
            spyOn(store.data,'getByKey').andCallThrough();

            col.initCanvases();

            expect(Ext.DomQuery.select).toHaveBeenCalled();
            expect(store.data.getByKey).toHaveBeenCalledWith('5');
            expect(col.initCanvas).toHaveBeenCalledWith(canvases[0].id, 150, 100, record.data.mol, record);
        })

        it('treestore issue #81', function() {
            // calls initCanvas on all canvas tags with class="x-chemdoodle-cols"
            // canvas id contains a record id
            var view = function() {};
            view.dom = function() { return 'dom'; };
            view.id = 'gridview-1234';
            var store = function() {};
            var record = { data:{ 'mol': caffeineMolFile } };
            store.getNodeById = function() { return record };
            var col = createDefault();
            col.grid = function() { };
            col.grid.getView = function() { return view };
            col.grid.getStore = function() { return store };
            var canvases = [{ id:'gridview-1234-5-2' }];

            spyOn(Ext.DomQuery,'select').andReturn(canvases);
            spyOn(col,'initCanvas');
            spyOn(store,'getNodeById').andCallThrough();

            col.initCanvases();

            expect(Ext.DomQuery.select).toHaveBeenCalled();
            expect(store.getNodeById).toHaveBeenCalledWith('5');
            expect(col.initCanvas).toHaveBeenCalledWith(canvases[0].id, 150, 100, record.data.mol, record);
        });
    });

    describe('initCanvas', function() {
        it('default', function() {
            var col = createDefault();
            var doodle = function() {};
            doodle.loadMolecule = function() {};
            spyOn(ChemDoodle, 'ViewerCanvas').andReturn(doodle);
            spyOn(ChemDoodle, 'readMOL').andReturn('caffiene');
            spyOn(doodle, 'loadMolecule');

            var record = { data:{ 'mol': caffeineMolFile } };
            col.initCanvas('gridview-1234-5-2', 150, 100, record.data.mol, record);

            expect(ChemDoodle.ViewerCanvas).toHaveBeenCalledWith('gridview-1234-5-2', 150, 100);
            expect(ChemDoodle.readMOL).toHaveBeenCalledWith(record.data.mol);
            expect(doodle.loadMolecule).toHaveBeenCalledWith('caffiene');
        });

        it('highlight', function() {
            var defspecs = new ChemDoodle.structures.VisualSpecifications();
            var hlspecs = new ChemDoodle.structures.VisualSpecifications();
            hlspecs.bonds_color = 'black';
            hlspecs.atoms_color = 'black';
            var col = Ext.create('Esc.chemdoodle.Column', {
                text: 'Molecule', dataIndex: 'mol', atomIndex:'atoms',
                canvasClass: 'x-chemdoodle-cols2',  width: 162,
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

            // mock ChemDoodle.ViewerCanvas
            var doodle = function() {};
            doodle.loadMolecule = function() {};
            doodle.specs = defspecs;
            spyOn(ChemDoodle, 'ViewerCanvas').andReturn(doodle);
            spyOn(doodle, 'loadMolecule');

            var record = { data:{ 'mol': caffeineMolFile, atoms: '1,2,3' } };
            col.initCanvas('gridview-1234-5-2', 150, 100, record.data.mol, record);

            expect(ChemDoodle.ViewerCanvas).toHaveBeenCalledWith('gridview-1234-5-2', 150, 100, true);
            expect(doodle.loadMolecule).toHaveBeenCalled();
            var mol = doodle.loadMolecule.mostRecentCall.args[0];
            expect(doodle.specs.atoms_color).toEqual('cyan');
            expect(doodle.specs.bonds_color).toEqual('cyan');
            expect(mol.atoms[0].specs).toBeUndefined();
            expect(mol.atoms[0].specs).toBeUndefined();
            expect(mol.atoms[1].specs).toEqual(hlspecs);
            expect(mol.atoms[1].specs).toEqual(hlspecs);
        });
    });
});
