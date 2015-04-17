describe('Molecules', function() {
  describe('store', function() {
    var store = null;
    var url = appRootBase + '/data/molecules.json';

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

    it('clearMzFilter, unfiltered', function() {
      spyOn(store, 'loadPage');

      store.clearMzFilter();

      expect(store.loadPage).not.toHaveBeenCalledWith(1);
      expect(store.getProxy().extraParams.mz).toBeUndefined();
    });

    it('clearMzFilter', function() {
      var mz = 1133;
      store.getProxy().extraParams.mz = mz;
      spyOn(store, 'loadPage');

      store.clearMzFilter();

      expect(store.loadPage).toHaveBeenCalledWith(1);
      expect(store.getProxy().extraParams.mz).toBeUndefined();
    });

    it('setMzFilter', function() {
      var mz = 1133;
      spyOn(store, 'loadPage');

      store.setMzFilter(mz);

      expect(store.loadPage).toHaveBeenCalledWith(1);
      expect(store.getProxy().extraParams.mz).toEqual(mz);
    });

    it('setPageSize', function() {
      spyOn(store, 'loadPage');
      expect(store.pageSize).toEqual(25); // default value

      store.setPageSize(10);

      expect(store.pageSize).toEqual(10);
      expect(store.loadPage).toHaveBeenCalledWith(1);
    });

    it('proxy exception', function() {
      spyOn(Ext.Error, 'handle').and.returnValue(true);

      var proxy = store.getProxy();
      proxy.fireEvent('exception', proxy, 'bla', 'foo');

      expect(Ext.Error.handle).toHaveBeenCalledWith({
        msg: 'Failed to load molecules from server',
        response: 'bla',
        operation: 'foo'
      });
    });

    describe('totalUnfilteredCount', function() {
      it('no rawdata', function() {
        expect(store.getTotalUnfilteredCount()).toEqual(0);
      });

      it('rawdata with totalUnfiltered', function() {
        store.getProxy().getReader().rawData = {
          totalUnfiltered: 123
        };
        expect(store.getTotalUnfilteredCount()).toEqual(123);
      });

      it('rawdata without totalUnfiltered', function() {
        store.getProxy().getReader().rawData = {};
        expect(store.getTotalUnfilteredCount()).toEqual(0);
      });
    });
  });

  describe('controller', function() {
    var store = null,
      ctrl = null;

    beforeEach(function(done) {
      if (!ctrl) {
        var app = Ext.create('Esc.magmaweb.resultsAppTest', {
          controllers: ['Molecules'],
          onBeforeLaunch: function() {
            // create the viewport components for this controller
            Ext.create('Esc.magmaweb.view.molecule.Panel');
          },
          launch: function() {
          }
        });
        ctrl = app.getController('Molecules');
      }

      if (!store) {
        store = ctrl.getStore('Molecules');
        // mock onLaunch of controller
        store.on('load', done, this, {single: true});
        store.load();
      } else {
        done();
      }
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
      expect(store.getProxy().url).toEqual(appRootBase + '/data/molecules.json');
    });

    it('scan filter', function() {
      expect(store.getProxy().extraParams).toEqual({});

      spyOn(store, 'loadPage');
      store.setScanFilter(1133);
      expect(store.getProxy().extraParams).toEqual({
        scanid: 1133
      });
      expect(store.loadPage).toHaveBeenCalledWith(1);

      store.removeScanFilter();
      expect(store.getProxy().extraParams).toEqual({});
      expect(store.loadPage).toHaveBeenCalledWith(1);
    });

    it('change page size', function() {
      spyOn(store, 'setPageSize');

      ctrl.onPageSizeChange({
        getValue: function() {
          return 123;
        }
      });

      expect(store.setPageSize).toHaveBeenCalledWith(123);
    });

    it('select molecule', function() {
      var record = store.getById(352);
      var rm = Ext.create('Ext.selection.RowModel');

      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').and.returnValue(false); // listeners dont hear any events
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

      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').and.returnValue(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      ctrl.onSelect(rm, record);
      ctrl.onDeselect(rm, record);

      expect(f.callback).toHaveBeenCalledWith('moleculedeselect', 352, jasmine.any(Object));
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('apply scan filter', function() {
      // mock list
      var list = {
        getSelectionModel: function() {
          return {
            hasSelection: function() {
              return false;
            }
          };
        }
      };
      spyOn(ctrl, 'getMoleculeList').and.returnValue(list);

      // mock store
      var mockedstore = {
        setScanFilter: function() {}
      };
      spyOn(ctrl, 'getMoleculesStore').and.returnValue(mockedstore);
      spyOn(mockedstore, 'setScanFilter');

      var scanid = 1133;
      ctrl.applyScanFilter(scanid);

      expect(mockedstore.setScanFilter).toHaveBeenCalledWith(scanid);
    });

    describe('clear scan filter', function() {
      var mockedstore, list, sm;
      beforeEach(function() {
        // mock selection model
        sm = jasmine.createSpyObj('selectionModel', ['clearSelections']);
        sm.hasSelection = function() {
          return false;
        };

        // mock list
        list = {
          getFragmentScoreColumn: function() {
            return scorecol;
          },
          getSelectionModel: function() {
            return sm;
          },
          clearMzFilter: function() {},
          hideFragmentScoreColumn: function() {}
        };

        spyOn(ctrl, 'getMoleculeList').and.returnValue(list);
        spyOn(list, 'clearMzFilter');

        // mock store
        mockedstore = {
          setScanFilter: function() {},
          removeScanFilter: function() {},
          sorters: new Ext.util.AbstractMixedCollection(false, function(item) {
            return item.id || item.property;
          }),
          clearMzFilter: jasmine.createSpy('clearMzFilter')
        };
        spyOn(ctrl, 'getMoleculesStore').and.returnValue(mockedstore);
        spyOn(mockedstore, 'setScanFilter');
        spyOn(mockedstore, 'removeScanFilter');
      });

      it('should remove scan filter from store', function() {
        var scanid = 1133;
        ctrl.applyScanFilter(scanid);

        ctrl.clearScanFilter();

        expect(mockedstore.removeScanFilter).toHaveBeenCalled();
      });

      it('should remove mz filter', function() {
        var scanid = 1133;
        ctrl.applyScanFilter(scanid);

        ctrl.clearScanFilter();

        expect(mockedstore.clearMzFilter).toHaveBeenCalled();
      });
    });

    it('clear filters', function() {
      // mock list
      var list = {
        clearFilters: function() {}
      };
      spyOn(ctrl, 'getMoleculeList').and.returnValue(list);
      spyOn(list, 'clearFilters');

      // mock application eventbus
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').and.returnValue(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      ctrl.clearFilters();

      expect(list.clearFilters).toHaveBeenCalled();
      expect(f.callback).toHaveBeenCalledWith('moleculenoselect');
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    describe('applyMzFilter', function() {
      var store = null,
        list = null,
        sm = null;
      beforeEach(function() {
        sm = {
          hasSelection: function() {
            return false;
          }
        };
        store = {
          isFilteredOnScan: function() {
            return true;
          },
          sorters: new Ext.util.AbstractMixedCollection(false, function(item) {
            return item.id || item.property;
          }),
          setMzFilter: jasmine.createSpy('setmzFilter')
        };
        spyOn(ctrl, 'getMoleculesStore').and.returnValue(store);
        list = jasmine.createSpyObj('list', ['showFragmentScoreColumn']);
        list.getSelectionModel = function() {
          return sm;
        };
        spyOn(ctrl, 'getMoleculeList').and.returnValue(list);
      });

      it('should not filter on mz of level > 1 scan', function() {
        ctrl.applyMzFilter(122.0373001, 2);

        expect(store.setMzFilter).not.toHaveBeenCalled();
      });

      it('should not filter on mz when no scan is selected', function() {
        store.isFilteredOnScan = function() {
          return false;
        };

        ctrl.applyMzFilter(122.0373001, 1);

        expect(store.setMzFilter).not.toHaveBeenCalled();
      });

      it('should tell list to filter on mz', function() {
        ctrl.applyMzFilter(122.0373001, 1);

        expect(store.setMzFilter).toHaveBeenCalledWith(122.0373001);
      });

      it('should show score column', function() {
        ctrl.applyMzFilter(122.0373001, 1);

        expect(list.showFragmentScoreColumn).toHaveBeenCalled();
      });

      it('should sort molecules on score, refscore and molid', function() {
        ctrl.applyMzFilter(122.0373001, 1);

        expect(store.sorters.indexOfKey('score')).toEqual(0);
        expect(store.sorters.getByKey('score').direction).toEqual('ASC');
        expect(store.sorters.indexOfKey('refscore')).toEqual(1);
        expect(store.sorters.getByKey('refscore').direction).toEqual('DESC');
        expect(store.sorters.indexOfKey('molid')).toEqual(2);
        expect(store.sorters.getByKey('molid').direction).toEqual('ASC');
      });
    });

    describe('clearMzFilter', function() {
      var list = null,
        sm = null,
        store = null;
      beforeEach(function() {
        sm = jasmine.createSpyObj('selectionModel', ['clearSelections']);
        sm.hasSelection = function() {
          return false;
        };
        list = jasmine.createSpyObj('list', ['clearMzFilter', 'hideFragmentScoreColumn']);
        list.getSelectionModel = function() {
          return sm;
        };
        spyOn(ctrl, 'getMoleculeList').and.returnValue(list);

        // mock store
        store = {
          setScanFilter: function() {},
          removeScanFilter: function() {},
          clearMzFilter: jasmine.createSpy('clearMzFilter'),
          sorters: new Ext.util.AbstractMixedCollection(false, function(item) {
            return item.id || item.property;
          })
        };
        spyOn(ctrl, 'getMoleculesStore').and.returnValue(store);
        spyOn(store, 'setScanFilter');
        spyOn(store, 'removeScanFilter');
      });

      it('should not filter on mz of level > 1 scan', function() {
        ctrl.clearMzFilter(122.0373001, 2);

        expect(store.clearMzFilter).not.toHaveBeenCalled();
      });

      it('should tell store to filter on mz', function() {
        ctrl.clearMzFilter(122.0373001, 1);

        expect(store.clearMzFilter).toHaveBeenCalledWith();
      });

      it('should hide score column', function() {
        ctrl.clearMzFilter(122.0373001, 1);

        expect(list.hideFragmentScoreColumn).toHaveBeenCalled();
      });

      it('should clear molecule selection when molecule is selected', function() {
        sm.hasSelection = function() { return true; };
        var callback = jasmine.createSpy('callback').and.returnValue(false);
        Ext.util.Observable.capture(ctrl.application, callback);

        ctrl.clearMzFilter(122.0373001, 1);

        expect(sm.clearSelections).toHaveBeenCalled();
        expect(callback).toHaveBeenCalledWith('moleculenoselect');

        Ext.util.Observable.releaseCapture(ctrl.application);
      });

      it('should not score filter and sort', function() {
        var scanid = 1133;
        ctrl.applyScanFilter(scanid);

        ctrl.clearScanFilter();

        expect(store.removeScanFilter).toHaveBeenCalled();
      });

      it('should clear score filter and sort', function() {
        store.sorters.add('score', [1, 2, 3]);
        store.filters = new Ext.util.MixedCollection();
        store.filters.add('score', [4, 5, 6]);

        ctrl.clearMzFilter(122.0373001, 1);

        expect(store.filters.containsKey('score')).toBeFalsy();
        expect(store.sorters.containsKey('score')).toBeFalsy();
      });
    });

    describe('onChromatrogramLoad', function() {
      var form;
      beforeEach(function() {
        form = {
          setDisabledAnnotateFieldset: function() {}
        };
        spyOn(form, 'setDisabledAnnotateFieldset');
        spyOn(ctrl, 'getMoleculeAddForm').and.returnValue(form);
      });

      it('initially', function() {
        expect(ctrl.hasMSData).toBeFalsy();
      });

      it('empty chromatogram', function() {
        chromatogram = {
          data: []
        };
        ctrl.onChromatrogramLoad(chromatogram);
        expect(ctrl.hasMSData).toBeFalsy();
        expect(form.setDisabledAnnotateFieldset).toHaveBeenCalledWith(true);
      });

      it('filled chromatogram', function() {
        chromatogram = {
          data: [1]
        };
        ctrl.onChromatrogramLoad(chromatogram);
        expect(ctrl.hasMSData).toBeTruthy();
        expect(form.setDisabledAnnotateFieldset).toHaveBeenCalledWith(false);
      });
    });

    it('load molecules', function() {
      spyOn(ctrl, 'metabolizable');
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').and.returnValue(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);
      var sm = { hasSelection: function() { return false; } };
      spyOn(ctrl, 'getSelectionModel').and.returnValue(sm);

      store.fireEvent('load', store, [], true);

      expect(f.callback).toHaveBeenCalledWith('moleculeload', jasmine.any(Object));
      expect(ctrl.metabolizable).toHaveBeenCalledWith(true);

      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('aborted load molecules', function() {
      spyOn(ctrl, 'metabolizable');

      ctrl.onLoad(store, [], false);

      expect(ctrl.metabolizable).not.toHaveBeenCalled();
    });

    it('load molecules, one molecule', function() {
      // mock list
      var sm = {
        hasSelection: function() {},
        select: function() {}
      };
      spyOn(ctrl, 'getSelectionModel').and.returnValue(sm);
      spyOn(sm, 'hasSelection').and.returnValue(false);
      spyOn(sm, 'select');

      spyOn(ctrl, 'metabolizable');
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').and.returnValue(false); // listeners dont hear any events

      Ext.util.Observable.capture(ctrl.application, f.callback);

      // fake loading a filtered list of molecules with only one molecule
      var record = store.getById(352);
      store.loadRecords([record]);
      store.fireEvent('load', store, [record], true);

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
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').and.returnValue(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);
      var sm = { hasSelection: function() { return false; } };
      spyOn(ctrl, 'getSelectionModel').and.returnValue(sm);

      // load zero
      store.loadRawData({
        rows: [],
        total: 0,
        totalUnfiltered: 0,
        scans: []
      });
      store.fireEvent('load', store, [], true);

      expect(f.callback).toHaveBeenCalledWith('moleculeload', jasmine.any(Object));
      expect(ctrl.metabolizable).toHaveBeenCalledWith(false);
      expect(ctrl.showAddStructuresForm).toHaveBeenCalledWith();

      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    describe('download molecules in csv', function() {
      it('default', function() {
        spyOn(window, 'open');
        spyOn(ctrl, 'getMoleculeList').and.returnValue({
          getFilterQuery: function() {
            return [];
          },
          getVisiblColumnIndices: function() {
            return ['name', 'score'];
          }
        });

        ctrl.download_csv();

        var url = Ext.urlAppend(appRootBase + '/data/molecules.csv', Ext.Object.toQueryString({
          page: 1,
          start: 0,
          limit: 10,
          sort: Ext.JSON.encode([{
            property: 'refscore',
            direction: 'DESC'
          }, {
            property: 'molid',
            direction: 'ASC'
          }]),
          cols: Ext.JSON.encode(['name', 'score'])
        }));
        expect(window.open).toHaveBeenCalledWith(url, 'moleculescsv');
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
        }]);
        spyOn(window, 'open');
        spyOn(ctrl, 'getMoleculeList').and.returnValue({
          getFilterQuery: function() {
            return {
              filter: filter
            };
          },
          getVisiblColumnIndices: function() {
            return ['name', 'score'];
          }
        });

        ctrl.download_csv();

        var url = Ext.urlAppend(appRootBase + '/data/molecules.csv', Ext.Object.toQueryString({
          scanid: 50,
          page: 1,
          start: 0,
          limit: 10,
          sort: Ext.JSON.encode([{
            property: 'refscore',
            direction: 'DESC'
          }, {
            property: 'molid',
            direction: 'ASC'
          }]),
          filter: filter,
          cols: Ext.JSON.encode(['name', 'score'])
        }));
        expect(window.open).toHaveBeenCalledWith(url, 'moleculescsv');
      });
    });

    it('showAddStructuresForm loads defaults and show it', function() {
      ctrl.hasMSData = false;
      var addform = {
        loadDefaults: function() {}
      };
      spyOn(addform, 'loadDefaults');
      spyOn(ctrl, 'getMoleculeAddForm').and.returnValue(addform);
      var panel = {
        setActiveItem: function() {}
      };
      spyOn(panel, 'setActiveItem');
      spyOn(ctrl, 'getMoleculePanel').and.returnValue(panel);

      ctrl.showAddStructuresForm();

      expect(panel.setActiveItem).toHaveBeenCalledWith(1);
      expect(addform.loadDefaults).toHaveBeenCalled();
    });

    it('showGrid', function() {
      var panel = {
        setActiveItem: function() {}
      };
      spyOn(panel, 'setActiveItem');
      spyOn(ctrl, 'getMoleculePanel').and.returnValue(panel);

      ctrl.showGrid();

      expect(panel.setActiveItem).toHaveBeenCalledWith(0);
    });

    it('addStructuresHandler', function() {
      var form = {
        isValid: function() {
          return true;
        },
        submit: function() {}
      };
      spyOn(form, 'submit');
      var panel = {
        getForm: function() {
          return form;
        }
      };
      spyOn(ctrl, 'getMoleculeAddForm').and.returnValue(panel);

      ctrl.addStructuresHandler();

      expect(form.submit).toHaveBeenCalledWith({
        url: appRootBase + '/rpc/' + ctrl.application.jobid + '/add_structures',
        submitEmptyText: false,
        waitMsg: jasmine.any(String),
        success: jasmine.any(Function),
        failure: jasmine.any(Function)
      });
    });

    it('metabolizeHandler', function() {
      var form = {
        isValid: function() {
          return true;
        },
        submit: function() {}
      };
      spyOn(form, 'submit');
      var wf = {
        getForm: function() {
          return form;
        }
      };
      ctrl.metabolizeForm = wf;

      ctrl.metabolizeHandler();

      expect(form.submit).toHaveBeenCalledWith({
        url: appRootBase + '/rpc/' + ctrl.application.jobid + '/metabolize',
        submitEmptyText: false,
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
        isValid: function() {
          return true;
        },
        submit: function() {}
      };
      spyOn(form, 'submit');
      var wf = {
        getForm: function() {
          return form;
        }
      };
      ctrl.metabolizeStructureForm = wf;

      ctrl.metabolizeOneHandler();

      expect(form.submit).toHaveBeenCalledWith({
        url: appRootBase + '/rpc/' + ctrl.application.jobid + '/metabolize_one',
        submitEmptyText: false,
        waitMsg: jasmine.any(String),
        success: jasmine.any(Function),
        failure: jasmine.any(Function)
      });
    });

    it('metabolizable', function() {
      var button = {
        setDisabled: function() {}
      };
      spyOn(button, 'setDisabled');
      spyOn(Ext, 'getCmp').and.returnValue(button);

      ctrl.metabolizable(true);

      expect(Ext.getCmp).toHaveBeenCalledWith('metabolizeaction');
      expect(button.setDisabled).toHaveBeenCalledWith(false);
    });

    it('showActionsMenu', function() {
      var event = {
        getXY: function() {
          return [5, 10];
        }
      };
      spyOn(ctrl.actionsMenu, 'showAt');

      ctrl.showActionsMenu('tool', event);

      expect(ctrl.actionsMenu.showAt).toHaveBeenCalledWith([5, 10]);
    });

    it('showDownloadMenu', function() {
      var event = {
        getXY: function() {
          return [5, 10];
        }
      };
      spyOn(ctrl.downloadMenu, 'showAt');

      ctrl.showDownloadMenu('tool', event);

      expect(ctrl.downloadMenu.showAt).toHaveBeenCalledWith([5, 10]);
    });

    it('has selection model', function() {
      var list = {
        getSelectionModel: function() {}
      };
      spyOn(ctrl, 'getMoleculeList').and.returnValue(list);
      spyOn(list, 'getSelectionModel');

      ctrl.getSelectionModel();

      expect(ctrl.getMoleculeList).toHaveBeenCalled();
      expect(list.getSelectionModel).toHaveBeenCalled();
    });

    describe('rememberSelectedMolecule', function() {
      var sm = null;
      var f = null;

      beforeEach(function() {
        sm = {
          hasSelection: function() {},
          getSelection: function() {
            return [{
              getId: function() {
                return 1234;
              }
            }];
          },
          select: function() {}
        };
        spyOn(ctrl, 'getSelectionModel').and.returnValue(sm);
      });

      it('has no selection and reselects nothing', function() {
        spyOn(sm, 'hasSelection').and.returnValue(false);

        ctrl.rememberSelectedMolecule();

        expect('molid' in store.getProxy().extraParams).toBeFalsy();
      });

      it('has selection and reselects', function() {
        spyOn(sm, 'hasSelection').and.returnValue(true);

        ctrl.rememberSelectedMolecule();

        expect(store.getProxy().extraParams.molid).toEqual(1234);
      });
    });

    it('onLaunch', function() {
      spyOn(ctrl, 'applyRole');
      var mocklist = {
        filters: {
          createFilters: function() {}
        },
        setPageSize: function() {},
        getSelectionModel: function() {}
      };
      spyOn(ctrl, 'getMoleculeList').and.returnValue(mocklist);

      ctrl.onLaunch();

      expect(ctrl.applyRole).toHaveBeenCalledWith();
      expect(ctrl.getMoleculeList).toHaveBeenCalledWith();
    });

    it('cantrun', function() {
      ctrl.application.canRun = false;
      assignbut = jasmine.createSpyObj('abut', ['hideCommandsColumn']);
      spyOn(ctrl, 'getMoleculeList').and.returnValue(assignbut);

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
