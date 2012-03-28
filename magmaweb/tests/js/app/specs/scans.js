describe('Scans controller', function() {
  var ctrl = null;
  var mocked_chromatogram = null;

  beforeEach(function() {
    if (!ctrl) {
      ctrl = Application.getController('Scans');
    }
    mocked_chromatogram = {
      cutoff: null,
      setLoading: function() {},
      setData: function() {},
      setExtractedIonChromatogram: function() {},
      clearScanSelection: function() {},
      selectScan: function() {},
      hasData: function() {},
      setMarkers: function() {},
      redraw: function() {},
      resetScales: function() {},
    };
    spyOn(ctrl, 'getChromatogram').andReturn(mocked_chromatogram);
  });

  it('onLaunch', function() {
    spyOn(mocked_chromatogram, 'setLoading');
    spyOn(d3, 'json');
    ctrl.onLaunch();
    expect(ctrl.getChromatogram()).toBeDefined();
    expect(mocked_chromatogram.cutoff).toEqual(2000000);
    expect(mocked_chromatogram.setLoading).toHaveBeenCalledWith(true);
    expect(d3.json).toHaveBeenCalledWith('data/chromatogram.json', jasmine.any(Function));
  });

  it('loadChromatogramCallback', function() {
    spyOn(mocked_chromatogram, 'setLoading');
    spyOn(mocked_chromatogram, 'setData');
    spyOn(ctrl, 'resetScans');
    var f = { callback: function() {} };
    spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
    Ext.util.Observable.capture(ctrl.application, f.callback);

    var data = [1,2,3,4];
    ctrl.loadChromatogramCallback(data);

    expect(mocked_chromatogram.setLoading).toHaveBeenCalledWith(false);
    expect(mocked_chromatogram.setData).toHaveBeenCalledWith(data);
    expect(ctrl.resetScans).toHaveBeenCalled();

    expect(f.callback).toHaveBeenCalledWith('chromatogramload', jasmine.any(Object));
    Ext.util.Observable.releaseCapture(ctrl.application);
  });

  it('loadChromatogramCallback error', function() {
    var oldhandle = Ext.Error.handle;
    Ext.Error.handle = function(err) {
        return true;
    };
    spyOn(Ext.Error, 'handle').andCallThrough();

    ctrl.loadChromatogramCallback(null);

    expect(Ext.Error.handle).toHaveBeenCalledWith({
        msg: 'Failed to load chromatogram from server',
        sourceMethod : 'loadChromatogramCallback', sourceClass : 'Esc.magmaweb.controller.Scans'
    });
    Ext.Error.handle = oldhandle;
  });

  it('clearExtractedIonChromatogram', function() {
    spyOn(mocked_chromatogram, 'setExtractedIonChromatogram');
    ctrl.clearExtractedIonChromatogram();
    expect(mocked_chromatogram.setExtractedIonChromatogram).toHaveBeenCalledWith([]);
  });

  it('loadExtractedIonChromatogram', function() {
    spyOn(d3, 'json');
    spyOn(mocked_chromatogram, 'setLoading');

    ctrl.loadExtractedIonChromatogram(352);

    expect(mocked_chromatogram.setLoading).toHaveBeenCalledWith(true);
    expect(d3.json).toHaveBeenCalledWith(
      'data/extractedionchromatogram.352.json',
      jasmine.any(Function)
    );
  });

  it('loadExtractedIonChromatogramCallback', function() {
    spyOn(mocked_chromatogram, 'setLoading');
    spyOn(mocked_chromatogram, 'setExtractedIonChromatogram');
    spyOn(ctrl, 'setScans');

    var data = {
        scans: [ 1,2,3 ],
        chromatogram: [ 4, 5, 6 ]
    };
    ctrl.loadExtractedIonChromatogramCallback(data);

    expect(mocked_chromatogram.setLoading).toHaveBeenCalledWith(false);
    expect(mocked_chromatogram.setExtractedIonChromatogram).toHaveBeenCalledWith(data.chromatogram);
    expect(ctrl.setScans).toHaveBeenCalledWith(data.scans);
  });

  it('loadExtractedIonChromatogramCallback error', function() {
      var oldhandle = Ext.Error.handle;
      Ext.Error.handle = function(err) {
          return true;
      };
      spyOn(Ext.Error, 'handle').andCallThrough();

      ctrl.loadExtractedIonChromatogramCallback(null);

      expect(Ext.Error.handle).toHaveBeenCalledWith({
          msg: 'Failed to load extracted ion chromatogram from server',
          sourceMethod : 'loadExtractedIonChromatogramCallback', sourceClass : 'Esc.magmaweb.controller.Scans'
      });
      Ext.Error.handle = oldhandle;
  });

  it('searchScan', function() {
    spyOn(Ext.MessageBox, 'prompt');
    ctrl.searchScan();
    expect(Ext.MessageBox.prompt).toHaveBeenCalledWith(
      jasmine.any(String),
      jasmine.any(String),
      jasmine.any(Function)
    );
  });

  it('clearScanSelection', function() {
    var f = { callback: function() {} };
    spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
    Ext.util.Observable.capture(ctrl.application, f.callback);
    spyOn(mocked_chromatogram, 'clearScanSelection');

    ctrl.clearScanSelection();

    expect(mocked_chromatogram.clearScanSelection).toHaveBeenCalled();
    expect(f.callback).toHaveBeenCalledWith('noselectscan');
    Ext.util.Observable.releaseCapture(ctrl.application);
  });

  it('selectScan', function() {
    var f = { callback: function() {} };
    spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
    Ext.util.Observable.capture(ctrl.application, f.callback);
    spyOn(mocked_chromatogram, 'selectScan');

    ctrl.selectScan(1133);

    expect(mocked_chromatogram.selectScan).toHaveBeenCalledWith(1133);
    expect(f.callback).toHaveBeenCalledWith('selectscan', 1133);
    Ext.util.Observable.releaseCapture(ctrl.application);
  });

  describe('setScansOfMetabolites', function() {

    it('filled', function() {
        // mock metabolite store
        var data = { rawData: { scans: [1,2] }};
        var proxy = { getReader: function() { return data; }};
        var store = { getProxy: function() { return proxy; }, getTotalCount: function() { return 1 }};
        spyOn(ctrl, 'setScans');

        ctrl.setScansOfMetabolites(store);

        expect(ctrl.scans_of_metabolites).toEqual([1, 2]);
        expect(ctrl.setScans).toHaveBeenCalledWith([1, 2]);
        expect(ctrl.hasStructures).toBeTruthy();
    });

    it('empty', function() {
        // mock metabolite store
        var data = { rawData: { scans: [] }};
        var proxy = { getReader: function() { return data; }};
        var store = { getProxy: function() { return proxy; }, getTotalCount: function() { return 0 }};
        spyOn(ctrl, 'setScans');

        ctrl.setScansOfMetabolites(store);

        expect(ctrl.scans_of_metabolites).toEqual([]);
        expect(ctrl.setScans).toHaveBeenCalledWith([]);
        expect(ctrl.hasStructures).toBeFalsy();
    });
  });

  it('resetScans', function() {
    spyOn(ctrl, 'setScans');
    ctrl.scans_of_metabolites = [1, 2];
    ctrl.resetScans();
    expect(ctrl.setScans).toHaveBeenCalledWith([1, 2]);
  });

  describe('setScans', function() {
    it('no chromatogram', function() {
      spyOn(mocked_chromatogram, 'hasData').andReturn(false);
      ctrl.setScans([]);
      expect(mocked_chromatogram.hasData).toHaveBeenCalled();
    });

    it('no scans', function() {
      spyOn(mocked_chromatogram, 'hasData').andReturn(true);

      var f = { callback: function() {} };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      ctrl.setScans([]);

      expect(f.callback).toHaveBeenCalledWith('noscansfound');
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('one scan, not prev selected', function() {
      spyOn(mocked_chromatogram, 'hasData').andReturn(true);
      spyOn(mocked_chromatogram, 'setMarkers');
      spyOn(mocked_chromatogram, 'selectScan');
      spyOn(ctrl, 'selectScan');
      var f = { callback: function() {} };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      var scans = [{
         "id" : 1133,
         "rt" : 1624.99
      }];
      ctrl.setScans(scans);

      expect(ctrl.selectScan).toHaveBeenCalledWith(1133);
      expect(mocked_chromatogram.setMarkers).toHaveBeenCalledWith(scans);
      expect(mocked_chromatogram.selectScan).not.toHaveBeenCalled();
      expect(f.callback).not.toHaveBeenCalledWith('noscansfound');
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('one scan, prev selected', function() {
      spyOn(mocked_chromatogram, 'hasData').andReturn(true);
      spyOn(mocked_chromatogram, 'setMarkers');
      spyOn(mocked_chromatogram, 'selectScan');
      spyOn(ctrl, 'selectScan');
      var f = { callback: function() {} };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      mocked_chromatogram.selectedScan = 1133;
      var scans = [{
         "id" : 1133,
         "rt" : 1624.99
      }];
      ctrl.setScans(scans);

      expect(ctrl.selectScan).not.toHaveBeenCalledWith(1133);
      expect(mocked_chromatogram.setMarkers).toHaveBeenCalledWith(scans);
      expect(mocked_chromatogram.selectScan).toHaveBeenCalledWith(1133);
      expect(f.callback).not.toHaveBeenCalledWith('noscansfound');
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('selected scan', function() {
      spyOn(mocked_chromatogram, 'hasData').andReturn(true);
      spyOn(mocked_chromatogram, 'setMarkers');
      spyOn(mocked_chromatogram, 'selectScan');
      spyOn(ctrl, 'selectScan');
      var f = { callback: function() {} };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      mocked_chromatogram.selectedScan = 1133;
      var scans = [{
        "id" : 1152,
        "rt" : 1652.25
      },{
       "id" : 1133,
       "rt" : 1624.99
      }];
      ctrl.setScans(scans);

      expect(ctrl.selectScan).not.toHaveBeenCalledWith(1133);
      expect(mocked_chromatogram.setMarkers).toHaveBeenCalledWith(scans);
      expect(mocked_chromatogram.selectScan).toHaveBeenCalledWith(1133);
      expect(f.callback).not.toHaveBeenCalledWith('noscansfound');
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('no selected scan', function() {
      spyOn(mocked_chromatogram, 'hasData').andReturn(true);
      spyOn(mocked_chromatogram, 'setMarkers');
      spyOn(mocked_chromatogram, 'selectScan');
      spyOn(ctrl, 'selectScan');
      var f = { callback: function() {} };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      var scans = [{
        "id" : 1152,
        "rt" : 1652.25
      },{
         "id" : 1133,
         "rt" : 1624.99
      }];
      ctrl.setScans(scans);

      expect(ctrl.selectScan).not.toHaveBeenCalledWith(1133);
      expect(mocked_chromatogram.setMarkers).toHaveBeenCalledWith(scans);
      expect(mocked_chromatogram.selectScan).not.toHaveBeenCalled();
      expect(f.callback).not.toHaveBeenCalledWith('noscansfound');
      Ext.util.Observable.releaseCapture(ctrl.application);
    });
  });

  it('center', function() {
      spyOn(mocked_chromatogram, 'resetScales');
      var tool = {
        up: function() {
            return {
              down: function() { return mocked_chromatogram; }
            };
        }
      };
      ctrl.center(tool);

      expect(mocked_chromatogram.resetScales).toHaveBeenCalled();
  });
});