describe('Metabolites', function() {
  describe('store', function() {
    var store = null;
    var url = 'data/metabolites.json';

    beforeEach(function() {
      if (!store) {
        store = Ext.create('Esc.magmaweb.store.Metabolites');
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
        var oldhandle = Ext.Error.handle;
        Ext.Error.handle = function(err) {
            return true;
        };
        spyOn(Ext.Error, 'handle').andCallThrough();

        var proxy = store.getProxy();
        proxy.fireEvent('exception', proxy, 'bla', 'foo');

        expect(Ext.Error.handle).toHaveBeenCalledWith({
            msg: 'Failed to load metabolites from server',
            response: 'bla',
            operation: 'foo'
        });
        Ext.Error.handle = oldhandle;
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
            ctrl = Application.getController('Metabolites');
        }

        if (!store) {
            store = ctrl.getStore('Metabolites');
            store.load(); // mock onLaunch of controller
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

     it('should have metabolites', function() {
        expect(store.getCount()).toBeGreaterThan(1);
     });

     it('controller configured store', function() {
        expect(store.getProxy().url).toEqual('data/metabolites.json');
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

     it('select metabolite', function() {
         var record = store.getById(352);
         var rm = Ext.create('Ext.selection.RowModel');

         var f = { callback: function() {} };
         spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
         Ext.util.Observable.capture(ctrl.application, f.callback);

         ctrl.onSelect(rm, record);
         //ctrl.getMetaboliteList().getSelectionModel().select([record]);

         expect(f.callback).toHaveBeenCalledWith('metaboliteselect', 352, jasmine.any(Object));
         Ext.util.Observable.releaseCapture(ctrl.application);
     });

     it('before select metabolite with scans', function() {
       var record = store.getById(352);
       var rm = Ext.create('Ext.selection.RowModel');
       expect(record.data.nr_scans).toBeGreaterThan(0);
       expect(ctrl.beforeSelect(rm, record)).toBeTruthy();
     });

     it('before select metabolite without scans', function() {
       var record = store.getById(78);
       var rm = Ext.create('Ext.selection.RowModel');
       expect(record.data.nr_scans).toEqual(0);
       expect(ctrl.beforeSelect(rm, record)).toBeFalsy();
     });

     it('deselect metabolite', function() {
       var record = store.getById(352);
       var rm = Ext.create('Ext.selection.RowModel');

       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
       Ext.util.Observable.capture(ctrl.application, f.callback);

       ctrl.onSelect(rm, record);
       ctrl.onDeselect(rm, record);

       expect(f.callback).toHaveBeenCalledWith('metabolitedeselect', 352, jasmine.any(Object));
       Ext.util.Observable.releaseCapture(ctrl.application);
     });

     it('apply scan filter', function() {
       // mock list
       var sm = { hasSelection: function() {} };
       var list = {
         getSelectionModel: function() { return sm; },
         showFragmentScoreColumn: function() {}
       };
       spyOn(ctrl, 'getMetaboliteList').andReturn(list);
       spyOn(list, 'showFragmentScoreColumn');
       spyOn(sm, 'hasSelection').andReturn(false);

       // mock store
       var mockedstore = { setScanFilter: function() {} };
       spyOn(ctrl, 'getMetabolitesStore').andReturn(mockedstore);
       spyOn(mockedstore, 'setScanFilter');

       var scanid = 1133;
       ctrl.applyScanFilter(scanid);

       expect(list.showFragmentScoreColumn).toHaveBeenCalled();
       expect(mockedstore.setScanFilter).toHaveBeenCalledWith(scanid);
       expect(sm.hasSelection).toHaveBeenCalled();
     });

     describe('clear scan filter', function() {
       var mockedstore, list, sm;
       beforeEach(function() {
           // mock list
           sm = { hasSelection: function() {}, deselectAll: function() {} };
           list = {
               getSelectionModel: function() { return sm; },
               getFragmentScoreColumn: function() { return scorecol; },
               hideFragmentScoreColumn: function() { },
               showFragmentScoreColumn: function() { }
           };

           spyOn(list, 'hideFragmentScoreColumn');
           spyOn(ctrl, 'getMetaboliteList').andReturn(list);
           spyOn(sm, 'hasSelection').andReturn(false);

           // mock store
           mockedstore = {
               setScanFilter: function() {},
               removeScanFilter: function() {}
           };
           spyOn(ctrl, 'getMetabolitesStore').andReturn(mockedstore);
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
           mockedstore.sorters = new Ext.util.MixedCollection();
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
       spyOn(ctrl, 'getMetaboliteList').andReturn(list);
       spyOn(list, 'clearFilters');

       // mock application eventbus
       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
       Ext.util.Observable.capture(ctrl.application, f.callback);

       ctrl.clearFilters();

       expect(list.clearFilters).toHaveBeenCalled();
       expect(f.callback).toHaveBeenCalledWith('metabolitenoselect');
       Ext.util.Observable.releaseCapture(ctrl.application);
     });

     describe('onChromatrogramLoad', function() {
        var form;
        beforeEach(function() {
            form = {
                setDisabledAnnotateFieldset: function() {}
            };
            spyOn(form, 'setDisabledAnnotateFieldset');
            spyOn(ctrl, 'getMetaboliteAddForm').andReturn(form);
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

     it('load metabolites', function() {
       spyOn(ctrl, 'metabolizable');
       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
       Ext.util.Observable.capture(ctrl.application, f.callback);

       store.fireEvent('load', store);

       expect(f.callback).toHaveBeenCalledWith('metaboliteload', jasmine.any(Object));
       expect(ctrl.metabolizable).toHaveBeenCalledWith(true);

       Ext.util.Observable.releaseCapture(ctrl.application);
     });

     it('load metabolites, one metabolite', function() {
       // mock list
       var sm = { hasSelection: function() {}, select: function() {} };
       var list = { getSelectionModel: function() { return sm; } };
       spyOn(ctrl, 'getMetaboliteList').andReturn(list);
       spyOn(sm, 'hasSelection').andReturn(false);
       spyOn(sm, 'select');

       spyOn(ctrl, 'metabolizable');
       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events

       Ext.util.Observable.capture(ctrl.application, f.callback);

       // fake loading a filtered list of metabolites with only one metabolite
       var record = store.getById(352);
       store.loadRecords([record]);
       store.fireEvent('load', store);

       expect(store.getCount()).toEqual(1);
       expect(sm.hasSelection).toHaveBeenCalled();
       expect(sm.select).toHaveBeenCalledWith(0);
       expect(f.callback).toHaveBeenCalledWith('metaboliteload', jasmine.any(Object));
       expect(ctrl.metabolizable).toHaveBeenCalledWith(true);
       Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('load metabolites, zero metabolites', function() {
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

       expect(f.callback).toHaveBeenCalledWith('metaboliteload', jasmine.any(Object));
       expect(ctrl.metabolizable).toHaveBeenCalledWith(false);
       expect(ctrl.showAddStructuresForm).toHaveBeenCalledWith();

       Ext.util.Observable.releaseCapture(ctrl.application);
    });

    describe('download metabolites', function() {
      it('default', function() {
          spyOn(window, 'open');

          // create fake filter
          spyOn(ctrl, 'getMetaboliteList').andReturn({
              getView: function() {
                  return {
                      getFeature: function() {
                          return {
                              getFilterData: function() {},
                              buildQuery: function() {
                                  return {};
                              }
                          }
                      }
                  }
              }
          })

          ctrl.download();
          var url = Ext.urlAppend('data/metabolites.csv', Ext.Object.toQueryString({
              page: 1,
              start: 0,
              limit: 10,
              sort: Ext.JSON.encode([{
                property: 'probability',
                direction: 'DESC'
              },{
                property: 'metid',
                direction: 'ASC'
              }])
          }));
          expect(window.open).toHaveBeenCalledWith(url ,'metabolites.csv');
      });

      it('filtered', function() {
          var proxy = store.getProxy();
          proxy.extraParams.scanid = 50;

          // create fake filter
          var filter = Ext.JSON.encode([{
              field: 'nr_scans',
              value: 1,
              type: 'numeric',
              comparison: 'gt'
          }])
          spyOn(ctrl, 'getMetaboliteList').andReturn({
              getView: function() {
                  return {
                      getFeature: function() {
                          return {
                              getFilterData: function() {},
                              buildQuery: function() {
                                  return { filter: filter };
                              }
                          }
                      }
                  }
              }
          })

          spyOn(window, 'open');
          ctrl.download();
          var url = Ext.urlAppend('data/metabolites.csv', Ext.Object.toQueryString({
              scanid: 50,
              page: 1,
              start: 0,
              limit: 10,
              sort: Ext.JSON.encode([{
                property: 'probability',
                direction: 'DESC'
              },{
                property: 'metid',
                direction: 'ASC'
              }]),
              filter: filter
          }));
          expect(window.open).toHaveBeenCalledWith(url ,'metabolites.csv');
      });
    });

    it('showAddStructuresForm loads defaults and show it', function() {
        ctrl.hasMSData = false;
        var addform = {
            loadDefaults: function() {}
        };
        spyOn(addform, 'loadDefaults');
        spyOn(ctrl, 'getMetaboliteAddForm').andReturn(addform);
        var panel = { setActiveItem: function() {} };
        spyOn(panel, 'setActiveItem');
        spyOn(ctrl, 'getMetabolitePanel').andReturn(panel);

        ctrl.showAddStructuresForm();

        expect(panel.setActiveItem).toHaveBeenCalledWith(1);
        expect(addform.loadDefaults).toHaveBeenCalled();
    });

    it('showGrid', function() {
        var panel = { setActiveItem: function() {} };
        spyOn(panel, 'setActiveItem');
        spyOn(ctrl, 'getMetabolitePanel').andReturn(panel);

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
        spyOn(ctrl, 'getMetaboliteAddForm').andReturn(panel);

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
            expect(ctrl.metabolizeForm.getForm().getValues()).toEqual({
                "metabolism_types": ["phase1", "phase2"],
                "n_reaction_steps": '1'
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
            setMetabolite: function() {},
            setDisabledAnnotateFieldset: function() {},
            show: function() {}
        };
        spyOn(form, 'setMetabolite');
        spyOn(form, 'setDisabledAnnotateFieldset');
        spyOn(form, 'show');
        ctrl.metabolizeStructureForm = form;

        ctrl.showMetabolizeStructureForm(1234);

        expect(form.setMetabolite).toHaveBeenCalledWith(1234);
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
  });
});