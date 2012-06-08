describe('MSpectras controller', function() {
  var ctrl = null;
  beforeEach(function() {
    if (!ctrl) {
      ctrl = Application.getController('MSpectras');
    }
  });

  it('config', function() {
    expect(ctrl.getMaxmslevel()).toEqual(3);
    expect(ctrl.getUrl()).toEqual('data/mspectra.{0}.json?mslevel={1}');
  });

  it('getMSpectra', function() {
    ctrl.mspectras[1] = 5;
    var mspectra = ctrl.getMSpectra(1);
    expect(mspectra).toEqual(5);
  });

  describe('loadMSpectra', function() {
    it('clearHigherSpectra==true', function() {
      var mslevel = 1, scanid = 1133, markers = [];
      var mspectra = Ext.create('Esc.d3.MSpectra');
      spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
      spyOn(mspectra, 'setLoading');
      spyOn(ctrl, 'clearMSpectraFrom');
      spyOn(d3, 'json');

      ctrl.loadMSpectra(mslevel, scanid, markers, true);

      expect(mspectra.setLoading).toHaveBeenCalledWith(true);
      expect(d3.json).toHaveBeenCalledWith(
        'data/mspectra.'+scanid+'.json?mslevel='+mslevel,
        jasmine.any(Function)
      );
      expect(ctrl.clearMSpectraFrom).toHaveBeenCalledWith(2);
      mspectra.destroy();
    });

    it('clearHigherSpectra==false', function() {
      var mslevel = 1, scanid = 1133, markers = [];
      var mspectra = Ext.create('Esc.d3.MSpectra');
      spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
      spyOn(mspectra, 'setLoading');
      spyOn(ctrl, 'clearMSpectraFrom');
      spyOn(d3, 'json');

      ctrl.loadMSpectra(mslevel, scanid, markers);

      expect(mspectra.setLoading).toHaveBeenCalledWith(true);
      expect(d3.json).toHaveBeenCalledWith(
        'data/mspectra.'+scanid+'.json?mslevel='+mslevel,
        jasmine.any(Function)
      );
      expect(ctrl.clearMSpectraFrom).not.toHaveBeenCalledWith(2);
      mspectra.destroy();
    });
  });

  describe('onLoadMSpectra', function() {
    it('notfound', function() {
      var mslevel = 1, scanid = 1133, markers = [];
      var mspectra = ctrl.getMSpectra(mslevel);
      var oldhandle = Ext.Error.handle;
      Ext.Error.handle = function(err) {
          return true;
      };
      spyOn(Ext.Error, 'handle').andCallThrough();

      ctrl.onLoadMSpectra(mslevel, scanid, markers, null);

      expect(Ext.Error.handle).toHaveBeenCalledWith({
          msg: 'Unable to find mspectra scan on level '+mslevel+' with id '+scanid,
          sourceMethod : 'onLoadMSpectra', sourceClass : 'Esc.magmaweb.controller.MSpectras'
      });
      Ext.Error.handle = oldhandle;
    });

    it('found', function() {
      var mslevel = 1, scanid = 1133, markers = [];
      var data = {
        cutoff: 200000,
        peaks: [{
          "mz" : 92.0256195068359,
          "intensity" : 969.546630859375
        },{
          "mz" : 93.1171875,
          "intensity" : 582.663879394531
        }]
      };
      var mspectra = Ext.create('Esc.d3.MSpectra');
      spyOn(mspectra, 'setLoading');
      spyOn(mspectra, 'setData');
      spyOn(mspectra, 'setMarkers');
      spyOn(mspectra, 'up').andCallFake(function() {
          return { down: function() {
              return { enable: function() { return true} };
          }};
      });
      spyOn(ctrl, 'getMSpectra').andReturn(mspectra);

      spyOn(Ext.MessageBox, 'show');
      var f = { callback: function() {} };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      ctrl.onLoadMSpectra(mslevel, scanid, markers, data);

      expect(mspectra.setLoading).toHaveBeenCalledWith(false);
      expect(mspectra.setData).toHaveBeenCalledWith(data.peaks);
      expect(mspectra.setMarkers).toHaveBeenCalledWith(markers);
      expect(mspectra.scanid).toEqual(scanid);
      expect(mspectra.cutoff).toEqual(data.cutoff);
      expect(mspectra.up).toHaveBeenCalled();
      expect(Ext.MessageBox.show).not.toHaveBeenCalled();
      expect(f.callback).toHaveBeenCalledWith('mspectraload', scanid, mslevel);
      Ext.util.Observable.releaseCapture(ctrl.application);
      mspectra.destroy();
    });
  });

  it('loadMSpectra2', function() {
    var scanid = 1134, markers = [1,2];
    spyOn(ctrl, 'loadMSpectra');
    ctrl.loadMSpectra2(scanid, markers);
    expect(ctrl.loadMSpectra).toHaveBeenCalledWith(2, scanid, markers, true);
  });

  it('loadMSpectra1', function() {
    var scanid = 1134;
    spyOn(ctrl, 'loadMSpectra');
    ctrl.loadMSpectra1(scanid);
    expect(ctrl.loadMSpectra).toHaveBeenCalledWith(1, scanid, [], true);
  });

  describe('loadMSpectraFromFragment', function() {
    it('already loaded', function() {
      var scanid = 1134, mslevel = 1;
      var child = { data: {scanid: scanid, mslevel: mslevel }};
      var frag = { firstChild: child, childNodes:[child] };
      spyOn(ctrl, 'getMSpectra').andReturn({ scanid: scanid});
      spyOn(ctrl, 'loadMSpectra');

      ctrl.loadMSpectraFromFragment(frag);

      expect(ctrl.getMSpectra).toHaveBeenCalledWith(mslevel);
      expect(ctrl.loadMSpectra).not.toHaveBeenCalled();
    });

    it('new scan', function() {
      var scanid = 1134, mslevel = 1;
      var child = { data: {scanid: scanid, mslevel: mslevel, mz: 123 }};
      var frag = { firstChild: child, childNodes:[child] };
      spyOn(ctrl, 'getMSpectra').andReturn({ scanid: null });
      spyOn(ctrl, 'loadMSpectra');

      ctrl.loadMSpectraFromFragment(frag);

      expect(ctrl.getMSpectra).toHaveBeenCalledWith(mslevel);
      expect(ctrl.loadMSpectra).toHaveBeenCalledWith(mslevel, scanid, [{mz:123}], true);
    });
  });

  describe('loadMSpectrasFromFragment', function() {
    it('lvl 1 fragment without children', function() {
      var scanid = 1134, mslevel = 1;
      var child = {
          data: {scanid: scanid, mslevel: mslevel, mz: 123 },
          hasChildNodes: function() { return false }
      };
      var frag = {
        firstChild: child,
        childNodes:[child],
        isRoot: function() {return true}
      };
      var mspectra = Ext.create('Esc.d3.MSpectra');
      spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
      spyOn(mspectra, 'setMarkers');
      spyOn(mspectra, 'selectPeak');
      spyOn(ctrl, 'loadMSpectra2');


      ctrl.loadMSpectrasFromFragment(frag, frag.childNodes);

      expect(mspectra.setMarkers).toHaveBeenCalledWith([{mz:123}]);
      expect(mspectra.selectPeak).toHaveBeenCalledWith(123);
      expect(ctrl.loadMSpectra2).not.toHaveBeenCalled();
      mspectra.destroy();
    });

    it('lvl 1 fragment with children', function() {
      var scanid = 1134, mslevel = 1;
      var child2 = { data: {scanid: 1135, mz: 56 }};
      var child = {
          data: {scanid: scanid, mslevel: mslevel, mz: 123 },
          hasChildNodes: function() { return true },
          firstChild: child2, childNodes: [child2]
      };
      var frag = {
        firstChild: child,
        childNodes:[child],
        isRoot: function() {return true}
      };
      var mspectra = Ext.create('Esc.d3.MSpectra');
      spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
      spyOn(mspectra, 'setMarkers');
      spyOn(mspectra, 'selectPeak');
      spyOn(ctrl, 'loadMSpectra2');

      ctrl.loadMSpectrasFromFragment(frag, frag.childNodes);

      expect(mspectra.setMarkers).toHaveBeenCalledWith([{mz:123}]);
      expect(mspectra.selectPeak).toHaveBeenCalledWith(123);
      expect(ctrl.loadMSpectra2).toHaveBeenCalledWith(
          1135,
          [{mz: 56}]
      );
      mspectra.destroy();
    });

    it('lvl 2 fragment', function() {
      var frag = {
        data: { mslevel:2, mz: 123 },
        isRoot: function() {return false}
      };
      var mspectra = Ext.create('Esc.d3.MSpectra');
      spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
      spyOn(mspectra, 'selectPeak');

      ctrl.loadMSpectrasFromFragment(frag, []);
      expect(mspectra.selectPeak).toHaveBeenCalledWith(123);
      mspectra.destroy();
    });
  });

  it('clearMSpectra', function() {
    var mslevel = 2;
    var mspectra = Ext.create('Esc.d3.MSpectra');
    spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
    spyOn(mspectra, 'setData');
    spyOn(mspectra, 'up').andCallFake(function() {
        return { down: function() {
            return { disable: function() { return true} };
        }};
    });
    var f = { callback: function() {} };
    spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
    Ext.util.Observable.capture(ctrl.application, f.callback);

    ctrl.clearMSpectra(mslevel);

    expect(mspectra.setData).toHaveBeenCalledWith([]);
    expect(mspectra.scanid).toEqual(-1);
    expect(f.callback).toHaveBeenCalledWith('mspectraclear', mslevel);
    expect(mspectra.up).toHaveBeenCalled();
    Ext.util.Observable.releaseCapture(ctrl.application);
    mspectra.destroy();
  });

  it('clearMSpectra1', function() {
    spyOn(ctrl, 'clearMSpectra');

    ctrl.clearMSpectra1();

    expect(ctrl.clearMSpectra).toHaveBeenCalledWith(1);
  });

  it('clearMSpectraFrom2', function() {
    var mspectra = Ext.create('Esc.d3.MSpectra');
    spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
    spyOn(ctrl, 'clearMSpectra');
    spyOn(mspectra, 'setMarkers');

    ctrl.clearMSpectraFrom2();

    expect(ctrl.getMSpectra(1).setMarkers).toHaveBeenCalledWith([]);
    expect(ctrl.clearMSpectra).not.toHaveBeenCalledWith(1);
    expect(ctrl.clearMSpectra).toHaveBeenCalledWith(2);
    expect(ctrl.clearMSpectra).toHaveBeenCalledWith(3);
    mspectra.destroy();
  });

  it('clearMSpectraFromFragment', function() {
    spyOn(ctrl, 'clearMSpectra');

    ctrl.clearMSpectraFromFragment({ data: { mslevel: 2}});

    expect(ctrl.clearMSpectra).not.toHaveBeenCalledWith(1);
    expect(ctrl.clearMSpectra).not.toHaveBeenCalledWith(2);
    expect(ctrl.clearMSpectra).toHaveBeenCalledWith(3);
  });

  describe('selectPeakFromFragment', function() {
    beforeEach(function() {
       function fakeMspectra() {
           return {
               selectPeak: function() {},
               clearPeakSelection: function() {},
           };
       }
       ctrl.mspectras[1] = fakeMspectra();
       ctrl.mspectras[2] = fakeMspectra();
       ctrl.mspectras[3] = fakeMspectra();
    });

    it('in lvl1', function() {
      spyOn(ctrl.mspectras[1], 'selectPeak');
      spyOn(ctrl.mspectras[2], 'clearPeakSelection');
      spyOn(ctrl.mspectras[3], 'clearPeakSelection');

      ctrl.selectPeakFromFragment({ data: {mslevel:1, mz: 123}});

      expect(ctrl.mspectras[1].selectPeak).toHaveBeenCalledWith(123);
      expect(ctrl.mspectras[2].clearPeakSelection).toHaveBeenCalled();
      expect(ctrl.mspectras[3].clearPeakSelection).toHaveBeenCalled();
    });

    it('in lvl2', function() {
      spyOn(ctrl.mspectras[1], 'selectPeak');
      spyOn(ctrl.mspectras[2], 'selectPeak');
      spyOn(ctrl.mspectras[3], 'clearPeakSelection');

      ctrl.selectPeakFromFragment({
        data: {mslevel:2, mz: 123},
        parentNode: { data: { mz: 456 }}
      });

      expect(ctrl.mspectras[2].selectPeak).toHaveBeenCalledWith(123);
      expect(ctrl.mspectras[1].selectPeak).toHaveBeenCalledWith(456);
      expect(ctrl.mspectras[3].clearPeakSelection).toHaveBeenCalled();
    });

    it('in lvl3', function() {
      spyOn(ctrl.mspectras[1], 'selectPeak');
      spyOn(ctrl.mspectras[2], 'selectPeak');
      spyOn(ctrl.mspectras[3], 'selectPeak');

      ctrl.selectPeakFromFragment({
        data: {mslevel:3, mz: 123},
        parentNode: {
          data: { mz: 456 },
          parentNode: { data: { mz: 789 } }
        }
      });

      expect(ctrl.mspectras[3].selectPeak).toHaveBeenCalledWith(123);
      expect(ctrl.mspectras[2].selectPeak).toHaveBeenCalledWith(456);
      expect(ctrl.mspectras[1].selectPeak).toHaveBeenCalledWith(789);
    });
  });

  it('deselectPeakFromFragment', function() {
    var mspectra = ctrl.getMSpectra(1);
    spyOn(mspectra, 'clearPeakSelection');

    ctrl.deselectPeakFromFragment({ data: {mslevel:1, mz: 123}});

    expect(mspectra.clearPeakSelection).toHaveBeenCalled();
  });

  it('center', function() {
      var mspectra = Ext.create('Esc.d3.MSpectra');
      spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
      spyOn(mspectra, 'resetScales');
      var tool = {
          up: function() {
              return {
                down: function() { return mspectra; }
              };
          }
      };

      ctrl.center(tool);

      expect(mspectra.resetScales).toHaveBeenCalled();
      mspectra.destroy();
  });
});