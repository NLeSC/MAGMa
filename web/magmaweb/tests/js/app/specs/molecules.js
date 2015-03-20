describe('Molecules', function() {
  describe('store', function() {
    var store = null;
    var url = 'data/molecules.json';

    beforeEach(function() {
      if (!store) {
        store = Ext.create('Esc.magmaweb.store.Molecules');
      }
    });

    it('create', function() {
      expect(store).toBeDefined();
      expect(store.getProxy().url).toBeUndefined();
      expect(store.pageSize).toEqual(25);
    });

    it('setUrl', function() {
      store.setUrl(url);
      expect(store.getProxy().url).toEqual(url);
    });

    it('removeScanFilter, unfiltered', function() {
      spyOn(store, 'loadPage');

      store.removeScanFilter();

      expect(store.loadPage).not.toHaveBeenCalledWith(1);
      expect(store.getProxy().extraParams.scanid).toBeUndefined();
    });

    it('removeScanFilter', function() {
      var scanid = 1133;
      store.getProxy().extraParams.scanid = scanid;
      spyOn(store, 'loadPage');

      store.removeScanFilter();

      expect(store.loadPage).toHaveBeenCalledWith(1);
      expect(store.getProxy().extraParams.scanid).toBeUndefined();
    });

    it('setScanFilter', function() {
      var scanid = 1133;
      spyOn(store, 'loadPage');

      store.setScanFilter(scanid);

      expect(store.loadPage).toHaveBeenCalledWith(1);
      expect(store.getProxy().extraParams.scanid).toEqual(scanid);
    });

    it('setPageSize', function() {
      spyOn(store, 'loadPage');
      expect(store.pageSize).toEqual(25); // default value

      store.setPageSize(10);

      expect(store.pageSize).toEqual(10);
      expect(store.loadPage).toHaveBeenCalledWith(1);
    });

    it('proxy exception', function() {
        spyOn(Ext.Error, 'handle').andReturn(true);

        var proxy = store.getProxy();
        proxy.fireEvent('exception', proxy, 'bla', 'foo');

        expect(Ext.Error.handle).toHaveBeenCalledWith({
            msg: 'Failed to load molecules from server',
            response: 'bla',
            operation: 'foo'
        });
    });

    describe('totalUnfilteredCount' , function() {
        it('no rawdata', function() {
            expect(store.getTotalUnfilteredCount()).toEqual(0);
        });

        it('rawdata with totalUnfiltered', function() {
            store.getProxy().getReader().rawData = {totalUnfiltered: 123};
            expect(store.getTotalUnfilteredCount()).toEqual(123);
        });

        it('rawdata without totalUnfiltered', function() {
            store.getProxy().getReader().rawData = {};
            expect(store.getTotalUnfilteredCount()).toEqual(0);
        });
    });
  });

  describe('controller', function() {
     var store = null, ctrl= null;

     beforeEach(function() {
        if (!ctrl) {
            ctrl = Application.getController('Molecules');
        }

        if (!store) {
            store = ctrl.getStore('Molecules');
            // disable reselecting selected row, tested in describe('reselect'
            store.removeListener('beforeLoad', ctrl.onBeforeLoad, ctrl);
            // mock onLaunch of controller
            store.load();
        }

        expect(store).toBeTruthy();

        waitsFor(
            function() { return !store.isLoading()},
            'load never completed',
            4000
        );
     });

     afterEach(function() {
        if (ctrl) {
            if ('metabolizeForm' in ctrl && 'destroy' in ctrl.metabolizeForm) {
                ctrl.metabolizeForm.destroy();
            }
            if ('metabolizeStructureForm' in ctrl && 'destroy' in ctrl.metabolizeStructureForm) {
                ctrl.metabolizeStructureForm.destroy();
            }
        }
     });

     it('should have molecules', function() {
        expect(store.getCount()).toBeGreaterThan(1);
     });

     it('controller configured store', function() {
        expect(store.getProxy().url).toEqual('data/molecules.json');
     });

     it('scan filter', function() {
       expect(store.getProxy().extraParams).toEqual({});

       spyOn(store, 'loadPage');
       store.setScanFilter(1133);
       expect(store.getProxy().extraParams).toEqual({ scanid: 1133});
       expect(store.loadPage).toHaveBeenCalledWith(1);

       store.removeScanFilter();
       expect(store.getProxy().extraParams).toEqual({});
       expect(store.loadPage).toHaveBeenCalledWith(1);
     });

     it('change page size', function() {
         spyOn(store, 'setPageSize');

         ctrl.onPageSizeChange({
             getValue: function() { return 123; }
         });

         expect(store.setPageSize).toHaveBeenCalledWith(123);
     });

     it('select molecule', function() {
         var record = store.getById(352);
         var rm = Ext.create('Ext.selection.RowModel');

         var f = { callback: function() {} };
         spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
         Ext.util.Observable.capture(ctrl.application, f.callback);

         ctrl.onSelect(rm, record);
         //ctrl.getMoleculeList().getSelectionModel().select([record]);

         expect(f.callback).toHaveBeenCalledWith('moleculeselect', 352, jasmine.any(Object));
         Ext.util.Observable.releaseCapture(ctrl.application);
     });

     it('before select molecule with scans', function() {
       var record = store.getById(352);
       var rm = Ext.create('Ext.selection.RowModel');
       expect(record.data.nhits).toBeGreaterThan(0);
       expect(ctrl.beforeSelect(rm, record)).toBeTruthy();
     });

     it('before select molecule without scans', function() {
       var record = store.getById(78);
       var rm = Ext.create('Ext.selection.RowModel');
       expect(record.data.nhits).toEqual(0);
       expect(ctrl.beforeSelect(rm, record)).toBeFalsy();
     });

     it('deselect molecule', function() {
       var record = store.getById(352);
       var rm = Ext.create('Ext.selection.RowModel');

       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
       Ext.util.Observable.capture(ctrl.application, f.callback);

       ctrl.onSelect(rm, record);
       ctrl.onDeselect(rm, record);

       expect(f.callback).toHaveBeenCalledWith('moleculedeselect', 352, jasmine.any(Object));
       Ext.util.Observable.releaseCapture(ctrl.application);
     });

     it('apply scan filter', function() {
       // mock list
       var list = {
         showFragmentScoreColumn: function() {}
       };
       spyOn(ctrl, 'getMoleculeList').andReturn(list);
       spyOn(list, 'showFragmentScoreColumn');

       // mock store
       var mockedstore = {
           setScanFilter: function() {},
           sorters: new Ext.util.AbstractMixedCollection(false, function(item) {
               return item.id || item.property;
           })
       };
       spyOn(ctrl, 'getMoleculesStore').andReturn(mockedstore);
       spyOn(mockedstore, 'setScanFilter');

       var scanid = 1133;
       ctrl.applyScanFilter(scanid);

       expect(list.showFragmentScoreColumn).toHaveBeenCalled();
       expect(mockedstore.setScanFilter).toHaveBeenCalledWith(scanid);
       expect(mockedstore.sorters.indexOfKey('score')).toEqual(0);
       expect(mockedstore.sorters.getByKey('score').direction).toEqual('ASC');
       expect(mockedstore.sorters.indexOfKey('refscore')).toEqual(1);
       expect(mockedstore.sorters.getByKey('refscore').direction).toEqual('DESC');
       expect(mockedstore.sorters.indexOfKey('molid')).toEqual(2);
       expect(mockedstore.sorters.getByKey('molid').direction).toEqual('ASC');
     });

     describe('clear scan filter', function() {
       var mockedstore, list, sm;
       beforeEach(function() {
           // mock list
           list = {
               getFragmentScoreColumn: function() { return scorecol; },
               hideFragmentScoreColumn: function() { },
               showFragmentScoreColumn: function() { }
           };

           spyOn(list, 'hideFragmentScoreColumn');
           spyOn(ctrl, 'getMoleculeList').andReturn(list);

           // mock store
           mockedstore = {
               setScanFilter: function() {},
               removeScanFilter: function() {},
               sorters: new Ext.util.AbstractMixedCollection(false, function(item) {
                   return item.id || item.property;
               })
           };
           spyOn(ctrl, 'getMoleculesStore').andReturn(mockedstore);
           spyOn(mockedstore, 'setScanFilter');
           spyOn(mockedstore, 'removeScanFilter');
       });

       it('without score filter and sort', function() {

           var scanid = 1133;
           ctrl.applyScanFilter(scanid);

           ctrl.clearScanFilter();

           expect(mockedstore.removeScanFilter).toHaveBeenCalled();
           expect(list.hideFragmentScoreColumn).toHaveBeenCalled();
       });

       it('with score filter and sort', function() {
           mockedstore.sorters.add('score', [1, 2, 3]);
           mockedstore.filters = new Ext.util.MixedCollection();
           mockedstore.filters.add('score', [4, 5, 6]);

           var scanid = 1133;
           ctrl.applyScanFilter(scanid);

           ctrl.clearScanFilter();

           expect(mockedstore.filters.containsKey('score')).toBeFalsy();
           expect(mockedstore.sorters.containsKey('score')).toBeFalsy();
       });
     });

     it('clear filters', function() {
       // mock list
       var list = { clearFilters: function() {} };
       spyOn(ctrl, 'getMoleculeList').andReturn(list);
       spyOn(list, 'clearFilters');

       // mock application eventbus
       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
       Ext.util.Observable.capture(ctrl.application, f.callback);

       ctrl.clearFilters();

       expect(list.clearFilters).toHaveBeenCalled();
       expect(f.callback).toHaveBeenCalledWith('moleculenoselect');
       Ext.util.Observable.releaseCapture(ctrl.application);
     });

     describe('onChromatrogramLoad', function() {
        var form;
        beforeEach(function() {
            form = {
                setDisabledAnnotateFieldset: function() {}
            };
            spyOn(form, 'setDisabledAnnotateFieldset');
            spyOn(ctrl, 'getMoleculeAddForm').andReturn(form);
        });

        it('initially', function() {
            expect(ctrl.hasMSData).toBeFalsy();
        });

        it('empty chromatogram', function() {
            chromatogram = { data: [] };
            ctrl.onChromatrogramLoad(chromatogram);
            expect(ctrl.hasMSData).toBeFalsy();
            expect(form.setDisabledAnnotateFieldset).toHaveBeenCalledWith(true);
        });

        it('filled chromatogram', function() {
            chromatogram = { data: [1] };
            ctrl.onChromatrogramLoad(chromatogram);
            expect(ctrl.hasMSData).toBeTruthy();
            expect(form.setDisabledAnnotateFieldset).toHaveBeenCalledWith(false);
        });
     });

     it('load molecules', function() {
       spyOn(ctrl, 'metabolizable');
       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
       Ext.util.Observable.capture(ctrl.application, f.callback);

       store.fireEvent('load', store);

       expect(f.callback).toHaveBeenCalledWith('moleculeload', jasmine.any(Object));
       expect(ctrl.metabolizable).toHaveBeenCalledWith(true);

       Ext.util.Observable.releaseCapture(ctrl.application);
     });

     it('load molecules, one molecule', function() {
       // mock list
       var sm = { hasSelection: function() {}, select: function() {} };
       spyOn(ctrl, 'getSelectionModel').andReturn(sm);
       spyOn(sm, 'hasSelection').andReturn(false);
       spyOn(sm, 'select');

       spyOn(ctrl, 'metabolizable');
       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events

       Ext.util.Observable.capture(ctrl.application, f.callback);

       // fake loading a filtered list of molecules with only one molecule
       var record = store.getById(352);
       store.loadRecords([record]);
       store.fireEvent('load', store);

       expect(store.getCount()).toEqual(1);
       expect(sm.hasSelection).toHaveBeenCalled();
       expect(sm.select).toHaveBeenCalledWith(0);
       expect(f.callback).toHaveBeenCalledWith('moleculeload', jasmine.any(Object));
       expect(ctrl.metabolizable).toHaveBeenCalledWith(true);
       Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('load molecules, zero molecules', function() {
       spyOn(ctrl, 'metabolizable');
       spyOn(ctrl, 'showAddStructuresForm');
       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
       Ext.util.Observable.capture(ctrl.application, f.callback);

       // load zero
       store.loadRawData({
            rows: [],
            total: 0,
            totalUnfiltered: 0,
            scans: []
       });
       store.fireEvent('load', store);

       expect(f.callback).toHaveBeenCalledWith('moleculeload', jasmine.any(Object));
       expect(ctrl.metabolizable).toHaveBeenCalledWith(false);
       expect(ctrl.showAddStructuresForm).toHaveBeenCalledWith();

       Ext.util.Observable.releaseCapture(ctrl.application);
    });

    describe('download molecules in csv', function() {
      it('default', function() {
          spyOn(window, 'open');
          spyOn(ctrl, 'getMoleculeList').andReturn({
             getFilterQuery: function() { return []},
             getVisiblColumnIndices: function() { return ['name', 'score']}
          })

          ctrl.download_csv();

          var url = Ext.urlAppend('data/molecules.csv', Ext.Object.toQueryString({
              page: 1,
              start: 0,
              limit: 10,
              sort: Ext.JSON.encode([{
                property: 'refscore',
                direction: 'DESC'
              },{
                property: 'molid',
                direction: 'ASC'
              }]),
              cols: Ext.JSON.encode(['name', 'score'])
          }));
          expect(window.open).toHaveBeenCalledWith(url ,'moleculescsv');
      });

      it('filtered', function() {
      	  // select scan
          var proxy = store.getProxy();
          proxy.extraParams.scanid = 50;

          // apply filter
          filter = Ext.JSON.encode([{
              field: 'nhits',
              value: 1,
              type: 'numeric',
              comparison: 'gt'
          }])
          spyOn(window, 'open');
          spyOn(ctrl, 'getMoleculeList').andReturn({
             getFilterQuery: function() { return {filter: filter}},
             getVisiblColumnIndices: function() { return ['name', 'score']}
          })

          ctrl.download_csv();

          var url = Ext.urlAppend('data/molecules.csv', Ext.Object.toQueryString({
              scanid: 50,
              page: 1,
              start: 0,
              limit: 10,
              sort: Ext.JSON.encode([{
                property: 'refscore',
                direction: 'DESC'
              },{
                property: 'molid',
                direction: 'ASC'
              }]),
              filter: filter,
              cols: Ext.JSON.encode(['name', 'score'])
          }));
          expect(window.open).toHaveBeenCalledWith(url ,'moleculescsv');
      });
    });

    it('showAddStructuresForm loads defaults and show it', function() {
        ctrl.hasMSData = false;
        var addform = {
            loadDefaults: function() {}
        };
        spyOn(addform, 'loadDefaults');
        spyOn(ctrl, 'getMoleculeAddForm').andReturn(addform);
        var panel = { setActiveItem: function() {} };
        spyOn(panel, 'setActiveItem');
        spyOn(ctrl, 'getMoleculePanel').andReturn(panel);

        ctrl.showAddStructuresForm();

        expect(panel.setActiveItem).toHaveBeenCalledWith(1);
        expect(addform.loadDefaults).toHaveBeenCalled();
    });

    it('showGrid', function() {
        var panel = { setActiveItem: function() {} };
        spyOn(panel, 'setActiveItem');
        spyOn(ctrl, 'getMoleculePanel').andReturn(panel);

        ctrl.showGrid();

        expect(panel.setActiveItem).toHaveBeenCalledWith(0);
    });

    it('addStructuresHandler', function() {
        var form = {
            isValid: function() { return true; },
            submit: function() {}
        };
        spyOn(form, 'submit');
        var panel = { getForm: function() { return form; } };
        spyOn(ctrl, 'getMoleculeAddForm').andReturn(panel);

        ctrl.addStructuresHandler();

        expect(form.submit).toHaveBeenCalledWith({
            url: '/rpc/'+Application.jobid+'/add_structures',
            submitEmptyText : false,
            waitMsg: jasmine.any(String),
            success: jasmine.any(Function),
            failure: jasmine.any(Function)
        });
    });

    it('showMetabolizeForm', function() {
        ctrl.showMetabolizeForm();

        waitsFor(
            function() { return !ctrl.metabolizeForm.loading;},
            'Form defaults never loaded',
            1000
        );

        expect(ctrl.metabolizeForm.isVisible()).toBeTruthy();
        ctrl.metabolizeForm.hide();

        runs(function() {
            // scenario gets serialized to json string
            expect(ctrl.metabolizeForm.getForm().getValues()).toEqual({
                "metabolize": "on",
                "scenario": '[{"type":"phase1","steps":2},{"type":"phase2","steps":1}]'
            });
        });
    });

    it('metabolizeHandler', function() {
        var form = {
            isValid: function() { return true; },
            submit: function() {}
        };
        spyOn(form, 'submit');
        var wf = { getForm: function() { return form }};
        ctrl.metabolizeForm = wf;

        ctrl.metabolizeHandler();

        expect(form.submit).toHaveBeenCalledWith({
            url: '/rpc/'+Application.jobid+'/metabolize',
            submitEmptyText : false,
            waitMsg: jasmine.any(String),
            success: jasmine.any(Function),
            failure: jasmine.any(Function)
        });
    });

    it('showMetabolizeStructureForm', function() {
        ctrl.hasMSData = false;
        var form = {
            setMolecule: function() {},
            setDisabledAnnotateFieldset: function() {},
            show: function() {}
        };
        spyOn(form, 'setMolecule');
        spyOn(form, 'setDisabledAnnotateFieldset');
        spyOn(form, 'show');
        ctrl.metabolizeStructureForm = form;

        ctrl.showMetabolizeStructureForm(1234);

        expect(form.setMolecule).toHaveBeenCalledWith(1234);
        expect(form.setDisabledAnnotateFieldset).toHaveBeenCalledWith(true);
        expect(form.show).toHaveBeenCalledWith();
    });

    it('metabolizeOneHandler', function() {
        var form = {
            isValid: function() { return true; },
            submit: function() {}
        };
        spyOn(form, 'submit');
        var wf = { getForm: function() { return form }};
        ctrl.metabolizeStructureForm = wf;

        ctrl.metabolizeOneHandler();

        expect(form.submit).toHaveBeenCalledWith({
            url: '/rpc/'+Application.jobid+'/metabolize_one',
            submitEmptyText : false,
            waitMsg: jasmine.any(String),
            success: jasmine.any(Function),
            failure: jasmine.any(Function)
        });
    });

    it('metabolizable', function() {
        var button = { setDisabled: function() {}};
        spyOn(button, 'setDisabled');
        spyOn(Ext, 'getCmp').andReturn(button);

        ctrl.metabolizable(true);

        expect(Ext.getCmp).toHaveBeenCalledWith('metabolizeaction');
        expect(button.setDisabled).toHaveBeenCalledWith(false);
    });

    it('showActionsMenu', function() {
        var event = { getXY: function() { return [5,10]; }};
        spyOn(ctrl.actionsMenu, 'showAt');

        ctrl.showActionsMenu('tool', event);

        expect(ctrl.actionsMenu.showAt).toHaveBeenCalledWith([5,10]);
    });

    it('showDownloadMenu', function() {
        var event = { getXY: function() { return [5,10]; }};
        spyOn(ctrl.downloadMenu, 'showAt');

        ctrl.showDownloadMenu('tool', event);

        expect(ctrl.downloadMenu.showAt).toHaveBeenCalledWith([5,10]);
    });

    it('has selection model', function() {
        var list = { getSelectionModel: function() {} };
        spyOn(ctrl, 'getMoleculeList').andReturn(list);
        spyOn(list, 'getSelectionModel');

        ctrl.getSelectionModel();

        expect(ctrl.getMoleculeList).toHaveBeenCalled();
        expect(list.getSelectionModel).toHaveBeenCalled();
    });

    describe('reselect', function() {
        var sm = null;
        var f = null;

        beforeEach(function() {
            sm = {
                hasSelection: function() {},
                getSelection: function() {
                    return [{getId: function() {}}];
                },
                select: function() {}
            };
            spyOn(ctrl, 'getSelectionModel').andReturn(sm);
        });

        it('has no selection and reselects nothing', function() {
            spyOn(sm, 'hasSelection').andReturn(false);
            spyOn(store, 'on');

            ctrl.onBeforeLoad(store);

            expect(store.on).not.toHaveBeenCalled();
        });

        it('has selection and reselects', function() {
            spyOn(sm, 'hasSelection').andReturn(true);
            spyOn(store, 'on');

            ctrl.onBeforeLoad(store);

            expect(store.on).toHaveBeenCalled();
        });
    });

    it('onLaunch', function() {
    	spyOn(ctrl, 'applyRole');
    	var mocklist = {
    			filters: {
    				createFilters: function() {}
    			},
    			setPageSize: function() {}
    	};
     	spyOn(ctrl, 'getMoleculeList').andReturn(mocklist);

    	ctrl.onLaunch();

    	expect(ctrl.applyRole).toHaveBeenCalledWith();
    	expect(ctrl.getMoleculeList).toHaveBeenCalledWith();
    });

    it('cantrun', function() {
  	  ctrl.application.canRun = false;
      assignbut = jasmine.createSpyObj('abut', [ 'hideCommandsColumn']);
      spyOn(ctrl,'getMoleculeList').andReturn(assignbut);

  	  ctrl.applyRole();

  	  expect(assignbut.hideCommandsColumn).toHaveBeenCalledWith();
    });

    it('showHelp', function() {
       spyOn(ctrl.application, 'showHelp');

       ctrl.showHelp();

       expect(ctrl.application.showHelp).toHaveBeenCalledWith('molpanel');
    });
  });
});