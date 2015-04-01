describe('MSpectras controller', function() {
  var ctrl = null;
  beforeEach(function() {
    if (!ctrl) {
      ctrl = Application.getController('MSpectras');
    }
  });

  function captureEvents() {
    var f = { callback: function() {} };
    spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
    Ext.util.Observable.capture(ctrl.application, f.callback);
    return f;
  }

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
	var mslevel = null, scanid = null, markers = null;
    var mspectra;

	beforeEach(function() {
	      mslevel = 1, scanid = 1133, markers = [];
	      mspectra = Ext.create('Esc.d3.MSpectra');
	      spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
	      spyOn(mspectra, 'setLoading');
	      spyOn(ctrl, 'clearMSpectraFrom');
	      spyOn(d3, 'json');
	});

	afterEach(function() {
		mspectra.destroy();
	});

    it('clearHigherSpectra==true', function() {
      ctrl.loadMSpectra(mslevel, scanid, markers, true);

      expect(mspectra.setLoading).toHaveBeenCalledWith(true);
      expect(d3.json).toHaveBeenCalledWith(
        'data/mspectra.'+scanid+'.json?mslevel='+mslevel,
        jasmine.any(Function)
      );
      expect(ctrl.clearMSpectraFrom).toHaveBeenCalledWith(2);
    });

    it('clearHigherSpectra==false', function() {
      ctrl.loadMSpectra(mslevel, scanid, markers);

      expect(mspectra.setLoading).toHaveBeenCalledWith(true);
      expect(d3.json).toHaveBeenCalledWith(
        'data/mspectra.'+scanid+'.json?mslevel='+mslevel,
        jasmine.any(Function)
      );
      expect(ctrl.clearMSpectraFrom).not.toHaveBeenCalledWith(2);
    });

    it('should not load spectra when it is already loaded', function() {
    	mspectra.scanid = scanid;

    	ctrl.loadMSpectra(mslevel, scanid, markers);

    	expect(d3.json).not.toHaveBeenCalled();
    });
  });

  describe('onLoadMSpectra', function() {
    it('notfound', function() {
      var mslevel = 1, scanid = 1133, markers = [];
      var mspectra = ctrl.getMSpectra(mslevel);
      spyOn(Ext.Error, 'handle').andReturn(true);

      ctrl.onLoadMSpectra(mslevel, scanid, markers, null);

      expect(Ext.Error.handle).toHaveBeenCalledWith({
          msg: 'Unable to find mspectra scan on level '+mslevel+' with id '+scanid,
          sourceMethod : 'onLoadMSpectra', sourceClass : 'Esc.magmaweb.controller.MSpectras'
      });
    });

    describe('lvl 1 spectra', function() {
      var mslevel = null, scanid = null;
      var markers = null, data = null;
      var mspectra = null, callback = null;

      beforeEach(function() {
    	  mslevel = 1;
    	  scanid = 1133;
          markers = [];
          data = {
            cutoff: 200000,
            peaks: [{
              "mz" : 92.0256195068359,
              "intensity" : 969.546630859375
            },{
              "mz" : 93.1171875,
              "intensity" : 582.663879394531
            }],
            "fragments": [{"mz": 92.0256195068359}, {"mz" : 93.1171875}]
          };
          mspectra = Ext.create('Esc.d3.MSpectra');
          spyOn(mspectra, 'setLoading');
          spyOn(mspectra, 'setData');
          spyOn(mspectra, 'setMarkers');
          spyOn(mspectra, 'selectPeak');
          spyOn(mspectra, 'up').andCallFake(function() {
              return { down: function() {
                  return { enable: function() { return true} };
              }};
          });
          spyOn(ctrl, 'getMSpectra').andReturn(mspectra);

          spyOn(Ext.MessageBox, 'show');
          callback = jasmine.createSpy('callback').andReturn(false); // listeners dont hear any events
          Ext.util.Observable.capture(ctrl.application, callback);
      });

      afterEach(function() {
          Ext.util.Observable.releaseCapture(ctrl.application);
          mspectra.destroy();
      });

      it('should render spectra', function() {
          ctrl.onLoadMSpectra(mslevel, scanid, markers, data);

          expect(mspectra.setLoading).toHaveBeenCalledWith(false);
          expect(mspectra.setData).toHaveBeenCalledWith(data.peaks);
          var expected_markers = [{"mz": 92.0256195068359}, {"mz" : 93.1171875}];
          expect(mspectra.setMarkers).toHaveBeenCalledWith(expected_markers);
          expect(mspectra.selectPeak).not.toHaveBeenCalled();
          expect(mspectra.scanid).toEqual(scanid);
          expect(mspectra.cutoff).toEqual(data.cutoff);
          expect(mspectra.up).toHaveBeenCalled();
          expect(Ext.MessageBox.show).not.toHaveBeenCalled();
          expect(callback).toHaveBeenCalledWith('mspectraload', scanid, mslevel);
      });

      it('should select peak when scan has only one peak', function() {
    	  data.fragments = [{"mz": 92.0256195068359}];

    	  ctrl.onLoadMSpectra(mslevel, scanid, markers, data);

          expect(mspectra.selectPeak).toHaveBeenCalledWith(92.0256195068359);
          expect(callback).toHaveBeenCalledWith('peakselect', 92.0256195068359, mslevel, scanid);
      });

      it('should select peak when peak was previously selected', function() {
    	  ctrl.selectedlvl1peak = 92.0256195068359;

    	  ctrl.onLoadMSpectra(mslevel, scanid, markers, data);

          expect(mspectra.selectPeak).toHaveBeenCalledWith(92.0256195068359);
          expect(callback).toHaveBeenCalledWith('peakselect', 92.0256195068359, mslevel, scanid);
      });
    });

    describe('lvl > 1 spectra', function() {
      it('should render spectra', function() {
          var mslevel = 2, scanid = 1133;
          var markers = [{"mz": 92.0256195068359}];
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
          var expected_markers = [{"mz": 92.0256195068359}];
          expect(mspectra.setMarkers).toHaveBeenCalledWith(expected_markers);
          expect(mspectra.scanid).toEqual(scanid);
          expect(mspectra.cutoff).toEqual(data.cutoff);
          expect(mspectra.up).toHaveBeenCalled();
          expect(Ext.MessageBox.show).not.toHaveBeenCalled();
          expect(f.callback).toHaveBeenCalledWith('mspectraload', scanid, mslevel);
          Ext.util.Observable.releaseCapture(ctrl.application);
          mspectra.destroy();
      });
    });
  });

  it('loadMSpectra2', function() {
    var scanid = 1134, markers = [1,2];
    spyOn(ctrl, 'loadMSpectra');
    ctrl.loadMSpectra2(scanid, markers);
    expect(ctrl.loadMSpectra).toHaveBeenCalledWith(2, scanid, markers, true);
  });

  describe('loadMSpectra1', function() {
	  beforeEach(function() {
		  spyOn(ctrl, 'loadMSpectra');
	  });

	  it('loadMSpectra1', function() {
	    var scanid = 1134;
	    ctrl.loadMSpectra1(scanid);

	    expect(ctrl.loadMSpectra).toHaveBeenCalledWith(1, scanid, [], true);
	  });

	  it('should clear selected peak', function() {
		ctrl.selectedlvl1peak = 311.2313;
		var mspectra = jasmine.createSpyObj('mspectra', ['clearPeakSelection']);
		mspectra.scanid = 1234;
		spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
		var callback = jasmine.createSpy('callback');
		Ext.util.Observable.capture(ctrl.application, callback);

		var scanid = 1134;
	    ctrl.loadMSpectra1(scanid);

	    expect(mspectra.clearPeakSelection).toHaveBeenCalledWith();
	    expect(callback).toHaveBeenCalledWith('peakdeselect', 311.2313, 1, 1134);
	    Ext.util.Observable.releaseCapture(ctrl.application);
	  });
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
      spyOn(ctrl, 'loadMSpectra2');

      ctrl.loadMSpectrasFromFragment(frag, frag.childNodes);

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
      spyOn(ctrl, 'loadMSpectra2');

      ctrl.loadMSpectrasFromFragment(frag, frag.childNodes);

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

  describe('clearMSpectra', function() {
	var mslevel = null;
	var mspectra = null;
	var callback = null;
	beforeEach(function() {
	    mslevel = 2;
	    mspectra = Ext.create('Esc.d3.MSpectra');
	    mspectra.scanid = 1234;
	    spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
	    spyOn(mspectra, 'setData');
	    spyOn(mspectra, 'up').andCallFake(function() {
	        return { down: function() {
	            return { disable: function() { return true} };
	        }};
	    });
	    callback = jasmine.createSpy('callback');
	    Ext.util.Observable.capture(ctrl.application, callback);
	});

	afterEach(function() {
	    Ext.util.Observable.releaseCapture(ctrl.application);
	    mspectra.destroy();
	});

	it('should clear mspectra', function() {
	    ctrl.clearMSpectra(mslevel);

	    expect(mspectra.setData).toHaveBeenCalledWith([]);
	    expect(mspectra.scanid).toEqual(false);
	    expect(callback).toHaveBeenCalledWith('mspectraclear', mslevel);
	    expect(mspectra.up).toHaveBeenCalled();
	});

	it('should not clear already clear mspectra', function() {
		mspectra.scanid = false;

		ctrl.clearMSpectra(mslevel);

		expect(mspectra.setData).not.toHaveBeenCalled();
	});
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
      expect(ctrl.mspectras[1].selectPeak).not.toHaveBeenCalled();
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
      expect(ctrl.mspectras[1].selectPeak).not.toHaveBeenCalled();
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


  it('showHelp', function() {
     spyOn(ctrl.application, 'showHelp');

     ctrl.showHelp();

     expect(ctrl.application.showHelp).toHaveBeenCalledWith('scan');
  });

  describe('selectPeakOfMolecule', function() {
	 it('should select peak with mz of molecule and fire event', function() {
		 var mol = {id: 1234, data: { mz: 5678.90}};
		 var mspectra = Ext.create('Esc.d3.MSpectra');
	     spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
	     spyOn(mspectra, 'selectPeak');
	     mspectra.scanid = 2;
	     var capturer = captureEvents();

		 ctrl.selectPeakOfMolecule(mol.id, mol);

		 expect(mspectra.selectPeak).toHaveBeenCalledWith(5678.90);
		 expect(capturer.callback).toHaveBeenCalledWith('peakselect', 5678.90, 1, 2);
		 Ext.util.Observable.releaseCapture(ctrl.application);
		 expect(ctrl.selectedlvl1peak).toEqual(5678.90);
	 });

	 it('should do nothing if peak already selected', function() {
		 var mol = {id: 1234, data: { mz: 5678.90}};
		 var mspectra = Ext.create('Esc.d3.MSpectra');
	     spyOn(ctrl, 'getMSpectra').andReturn(mspectra);
	     spyOn(mspectra, 'selectPeak');
	     mspectra.selectedpeak = 5678.90;
	     var capturer = captureEvents();

		 ctrl.selectPeakOfMolecule(mol.id, mol);

		 expect(mspectra.selectPeak).not.toHaveBeenCalled();
		 expect(capturer.callback).not.toHaveBeenCalled();
		 Ext.util.Observable.releaseCapture(ctrl.application);
	 });
  });

  describe('deselectPeakOfMolecule', function() {
	 it('should unselect lvl1 peak', function() {
		 ctrl.selectedlvl1peak = 5678.90;

		 ctrl.deselectPeakOfMolecule();

		 expect(ctrl.selectedlvl1peak).toBeFalsy();
	 }) ;
  });
});