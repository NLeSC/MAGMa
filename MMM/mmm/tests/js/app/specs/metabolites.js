describe('Metabolites', function() {
  describe('store', function() {
    var store = null;
    var url = 'data/metabolites.json';

    beforeEach(function() {
      if (!store) {
        store = Ext.create('Esc.mmm.store.Metabolites');
      }
    });

    it('create', function() {
      expect(store).toBeDefined();
      expect(store.getProxy().url).toBeUndefined();
      expect(store.isLoaded).toBeFalsy();
    });

    it('setUrl', function() {
      store.setUrl(url);
      expect(store.getProxy().url).toEqual(url);
    });

    it('isLoaded true', function() {
      store.fireEvent('load', store);
      expect(store.isLoaded).toBeTruthy();
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
       var scorecol = { show: function() {} };
       var list = {
         getSelectionModel: function() { return sm; },
         getFragmentScoreColumn: function() { return scorecol; }
       };
       spyOn(ctrl, 'getMetaboliteList').andReturn(list);
       spyOn(scorecol, 'show');
       spyOn(sm, 'hasSelection').andReturn(false);

       // mock store
       var mockedstore = { setScanFilter: function() {} };
       spyOn(ctrl, 'getMetabolitesStore').andReturn(mockedstore);
       spyOn(mockedstore, 'setScanFilter');

       var scanid = 1133;
       ctrl.applyScanFilter(scanid);

       expect(scorecol.show).toHaveBeenCalled();
       expect(mockedstore.setScanFilter).toHaveBeenCalledWith(scanid);
       expect(sm.hasSelection).toHaveBeenCalled();
     });

     describe('clear scan filter', function() {
       var mockedstore, scorecol, sm;
       beforeEach(function() {
           // mock list
           sm = { hasSelection: function() {} };
           scorecol = { show: function() {}, hide: function() {} };
           var list = {
               getSelectionModel: function() { return sm; },
               getFragmentScoreColumn: function() { return scorecol; }
           };
           spyOn(ctrl, 'getMetaboliteList').andReturn(list);
           spyOn(sm, 'hasSelection').andReturn(false);

           // mock store
           mockedstore = {
               setScanFilter: function() {},
               removeScanFilter: function() {}
           };
           spyOn(ctrl, 'getMetabolitesStore').andReturn(mockedstore);
           spyOn(scorecol, 'hide');
           spyOn(mockedstore, 'setScanFilter');
           spyOn(mockedstore, 'removeScanFilter');
       });

       it('without score filter and sort', function() {

           var scanid = 1133;
           ctrl.applyScanFilter(scanid);

           ctrl.clearScanFilter();

           expect(scorecol.hide).toHaveBeenCalled();
           expect(mockedstore.removeScanFilter).toHaveBeenCalled();
           expect(sm.hasSelection).toHaveBeenCalled();
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

     it('load metabolites', function() {
       var f = { callback: function() {} };
       spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
       Ext.util.Observable.capture(ctrl.application, f.callback);

       store.fireEvent('load', store);

       expect(f.callback).toHaveBeenCalledWith('metaboliteload', jasmine.any(Object));

       Ext.util.Observable.releaseCapture(ctrl.application);
     });

     it('load metabolites, one metabolite', function() {
       // mock list
       var sm = { hasSelection: function() {}, select: function() {} };
       var list = { getSelectionModel: function() { return sm; } };
       spyOn(ctrl, 'getMetaboliteList').andReturn(list);
       spyOn(sm, 'hasSelection').andReturn(false);
       spyOn(sm, 'select');

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
       Ext.util.Observable.releaseCapture(ctrl.application);
    });

  });
});