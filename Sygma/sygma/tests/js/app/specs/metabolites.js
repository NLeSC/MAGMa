describe('Metabolites', function() {
   var store = null, ctrl= null;

   beforeEach(function() {
      if (!ctrl) {
          ctrl = Application.getController('Metabolites');
      }

      if (!store) {
          store = ctrl.getStore('Metabolites');
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
      console.log(ctrl.getMetaboliteList());
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
       var f = { callback: function(metid, metabolite) {} };
       spyOn(f, 'callback');
       ctrl.application.addListener('metaboliteselect', f.callback);

       ctrl.onSelect(rm, record);
       expect(f.callback).toHaveBeenCalled();
       expect(f.callback.mostRecentCall.args[0]).toEqual(352);

       ctrl.application.removeListener('metaboliteselect', f.callback);
   });
});