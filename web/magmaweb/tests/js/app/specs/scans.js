describe('Scans controller', function() {
  var ctrl = null;
  var mocked_chromatogram = null;
  var mocked_form_panel = null;
  var mocked_form = null;

  beforeEach(function() {
    if (!ctrl) {
      var app = Ext.create('Esc.magmaweb.resultsAppTest', {
        controllers: ['Scans'],
        onBeforeLaunch: function() {
          this.initQuickTips();
          this.initViewport();
          // prevent controllers onLaunch from fireing as it will do an ajax call
          // ajax call will be tested later
        },
        launch: function() {
          // for some reason the view port is not auto created so do it here
        }
      });
      ctrl = app.getController('Scans');
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
      setZoom: function() {}
    };
    spyOn(ctrl, 'getChromatogram').andReturn(mocked_chromatogram);
    mocked_form = {
      submit: function() {},
      isValid: function() {
        return true;
      },
      setValues: function() {},
      getValues: function() {},
      load: function() {}
    };
    mocked_form_panel = {
      setDisabledAnnotateFieldset: function() {},
      loadDefaults: function() {},
      getForm: function() {
        return mocked_form;
      }
    };
    spyOn(ctrl, 'getUploadForm').andReturn(mocked_form_panel);
  });

  it('onLaunch', function() {
    spyOn(mocked_chromatogram, 'setLoading');
    spyOn(d3, 'json');
    spyOn(ctrl, 'applyRole');

    ctrl.onLaunch();

    expect(ctrl.getChromatogram()).toBeDefined();
    expect(mocked_chromatogram.setLoading).toHaveBeenCalledWith(true);
    expect(d3.json).toHaveBeenCalledWith(appRootBase + '/data/chromatogram.json', jasmine.any(Function));
    expect(ctrl.applyRole).toHaveBeenCalledWith();
  });

  describe('loadChromatogramCallback', function() {
    it('with cutoff', function() {
      spyOn(mocked_chromatogram, 'setLoading');
      spyOn(mocked_chromatogram, 'setData');
      spyOn(ctrl, 'resetScans');
      var mocked_panel = jasmine.createSpyObj('panel', ['hide']);
      spyOn(ctrl, 'getChromatogramPanel').andReturn(mocked_panel);
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      var data = {
        scans: [1, 2, 3, 4],
        cutoff: 1000
      };
      ctrl.loadChromatogramCallback(data);

      expect(mocked_chromatogram.setLoading).toHaveBeenCalledWith(false);
      expect(mocked_chromatogram.cutoff).toEqual(1000);
      expect(mocked_chromatogram.setData).toHaveBeenCalledWith([1, 2, 3, 4]);
      expect(ctrl.resetScans).toHaveBeenCalled();
      expect(mocked_panel.hide).not.toHaveBeenCalledWith();
      expect(f.callback).toHaveBeenCalledWith('chromatogramload', jasmine.any(Object));
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('with default cutoff', function() {
      mocked_chromatogram.cutoff = 1000;
      spyOn(mocked_chromatogram, 'setLoading');
      spyOn(mocked_chromatogram, 'setData');
      spyOn(ctrl, 'resetScans');
      var mocked_panel = jasmine.createSpyObj('panel', ['hide']);
      spyOn(ctrl, 'getChromatogramPanel').andReturn(mocked_panel);
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      var data = {
        scans: [1, 2, 3, 4],
        cutoff: null
      };
      ctrl.loadChromatogramCallback(data);

      expect(mocked_chromatogram.setLoading).toHaveBeenCalledWith(false);
      expect(mocked_chromatogram.cutoff).toEqual(1000);
      expect(mocked_chromatogram.setData).toHaveBeenCalledWith([1, 2, 3, 4]);
      expect(ctrl.resetScans).toHaveBeenCalled();
      expect(mocked_panel.hide).not.toHaveBeenCalledWith();
      expect(f.callback).toHaveBeenCalledWith('chromatogramload', jasmine.any(Object));
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('server error', function() {
      spyOn(Ext.Error, 'handle').andReturn(true);

      ctrl.loadChromatogramCallback(null);

      expect(Ext.Error.handle).toHaveBeenCalledWith({
        msg: 'Failed to load chromatogram from server',
        sourceMethod: 'loadChromatogramCallback',
        sourceClass: 'Esc.magmaweb.controller.Scans'
      });
    });

    it('when has a single scan then collapse chromatogram', function() {
      spyOn(mocked_chromatogram, 'setLoading');
      spyOn(mocked_chromatogram, 'setData');
      spyOn(mocked_chromatogram, 'selectScan');
      spyOn(ctrl, 'resetScans');
      var mocked_panel = jasmine.createSpyObj('panel', ['collapse']);
      spyOn(ctrl, 'getChromatogramPanel').andReturn(mocked_panel);
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      var data = {
        scans: [{
          id: 1234
        }],
        cutoff: null
      };
      ctrl.loadChromatogramCallback(data);

      expect(mocked_panel.collapse).toHaveBeenCalledWith();
      expect(mocked_chromatogram.selectScan).toHaveBeenCalledWith(1234, false);
      expect(f.callback).toHaveBeenCalledWith('chromatogramload', jasmine.any(Object));
      expect(f.callback).toHaveBeenCalledWith('selectscan', 1234);
      Ext.util.Observable.releaseCapture(ctrl.application);
    });
  });

  describe('ExtractedIonChromatogram', function() {
    it('clear', function() {
      spyOn(mocked_chromatogram, 'setExtractedIonChromatogram');
      ctrl.clearExtractedIonChromatogram();
      expect(mocked_chromatogram.setExtractedIonChromatogram).toHaveBeenCalledWith([]);
    });

    it('load', function() {
      spyOn(d3, 'json');
      spyOn(mocked_chromatogram, 'setLoading');

      ctrl.loadExtractedIonChromatogram(352);

      expect(mocked_chromatogram.setLoading).toHaveBeenCalledWith(true);
      expect(d3.json).toHaveBeenCalledWith(
        appRootBase + '/data/extractedionchromatogram.352.json',
        jasmine.any(Function)
      );
    });

    it('load callback success', function() {
      spyOn(mocked_chromatogram, 'setLoading');
      spyOn(mocked_chromatogram, 'setExtractedIonChromatogram');
      spyOn(ctrl, 'setScans');

      var data = {
        scans: [1, 2, 3],
        chromatogram: [4, 5, 6]
      };
      ctrl.loadExtractedIonChromatogramCallback(data);

      expect(mocked_chromatogram.setLoading).toHaveBeenCalledWith(false);
      expect(mocked_chromatogram.setExtractedIonChromatogram).toHaveBeenCalledWith(data.chromatogram);
      expect(ctrl.setScans).toHaveBeenCalledWith(data.scans);
    });

    it('load callback failure', function() {
      spyOn(Ext.Error, 'handle').andReturn(true);

      ctrl.loadExtractedIonChromatogramCallback(null);

      expect(Ext.Error.handle).toHaveBeenCalledWith({
        msg: 'Failed to load extracted ion chromatogram from server',
        sourceMethod: 'loadExtractedIonChromatogramCallback',
        sourceClass: 'Esc.magmaweb.controller.Scans'
      });
    });

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
    var f = {
      callback: function() {}
    };
    spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
    Ext.util.Observable.capture(ctrl.application, f.callback);
    spyOn(mocked_chromatogram, 'clearScanSelection');

    ctrl.clearScanSelection();

    expect(mocked_chromatogram.clearScanSelection).toHaveBeenCalled();
    expect(f.callback).toHaveBeenCalledWith('noselectscan');
    Ext.util.Observable.releaseCapture(ctrl.application);
  });

  describe('selectScan', function() {
    var f = null;

    beforeEach(function() {
      f = {
        callback: function() {}
      };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);
      spyOn(mocked_chromatogram, 'selectScan');
    });

    afterEach(function() {
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('default', function() {
      ctrl.selectScan(1133);

      expect(mocked_chromatogram.selectScan).toHaveBeenCalledWith(1133, undefined);
      expect(f.callback).toHaveBeenCalledWith('selectscan', 1133);
    });

    it('silent', function() {
      ctrl.selectScan(1133, true);

      expect(mocked_chromatogram.selectScan).toHaveBeenCalledWith(1133, true);
      expect(f.callback).toHaveBeenCalledWith('selectscan', 1133);
    });
  });

  describe('setScansOfMolecules', function() {

    it('filled', function() {
      // mock molecule store
      var data = {
        rawData: {
          scans: [1, 2]
        }
      };
      var proxy = {
        getReader: function() {
          return data;
        }
      };
      var store = {
        getProxy: function() {
          return proxy;
        },
        getTotalCount: function() {
          return 1;
        }
      };
      spyOn(ctrl, 'setScans');

      ctrl.setScansOfMolecules(store);

      expect(ctrl.scans_of_molecules).toEqual([1, 2]);
      expect(ctrl.setScans).toHaveBeenCalledWith([1, 2]);
      expect(ctrl.hasStructures).toBeTruthy();
    });

    it('empty', function() {
      // mock molecule store
      var data = {
        rawData: {
          scans: []
        }
      };
      var proxy = {
        getReader: function() {
          return data;
        }
      };
      var store = {
        getProxy: function() {
          return proxy;
        },
        getTotalCount: function() {
          return 0;
        }
      };
      spyOn(ctrl, 'setScans');

      ctrl.setScansOfMolecules(store);

      expect(ctrl.scans_of_molecules).toEqual([]);
      expect(ctrl.setScans).toHaveBeenCalledWith([]);
      expect(ctrl.hasStructures).toBeFalsy();
    });
  });

  it('resetScans', function() {
    spyOn(ctrl, 'setScans');
    ctrl.scans_of_molecules = [1, 2];
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

      var f = {
        callback: function() {}
      };
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
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      var scans = [{
        "id": 1133,
        "rt": 1624.99
      }];
      ctrl.setScans(scans);

      expect(ctrl.selectScan).not.toHaveBeenCalledWith(1133);
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
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      mocked_chromatogram.selectedScan = 1133;
      var scans = [{
        "id": 1133,
        "rt": 1624.99
      }];
      ctrl.setScans(scans);

      expect(ctrl.selectScan).not.toHaveBeenCalledWith(1133);
      expect(mocked_chromatogram.setMarkers).toHaveBeenCalledWith(scans);
      expect(mocked_chromatogram.selectScan).toHaveBeenCalledWith(1133);
      expect(f.callback).not.toHaveBeenCalledWith('noscansfound');
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('one scan, prev other scan selected', function() {
      spyOn(mocked_chromatogram, 'hasData').andReturn(true);
      spyOn(mocked_chromatogram, 'setMarkers');
      spyOn(mocked_chromatogram, 'clearScanSelection');
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      mocked_chromatogram.selectedScan = 204;
      var scans = [{
        "id": 1133,
        "rt": 1624.99
      }, {
        "id": 1134,
        "rt": 1634.99
      }];
      ctrl.setScans(scans);

      expect(mocked_chromatogram.setMarkers).toHaveBeenCalledWith(scans);
      expect(mocked_chromatogram.clearScanSelection).not.toHaveBeenCalled();
      expect(f.callback).not.toHaveBeenCalledWith('noselectscan');
      Ext.util.Observable.releaseCapture(ctrl.application);
    });

    it('selected scan', function() {
      spyOn(mocked_chromatogram, 'hasData').andReturn(true);
      spyOn(mocked_chromatogram, 'setMarkers');
      spyOn(mocked_chromatogram, 'selectScan');
      spyOn(ctrl, 'selectScan');
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      mocked_chromatogram.selectedScan = 1133;
      var scans = [{
        "id": 1152,
        "rt": 1652.25
      }, {
        "id": 1133,
        "rt": 1624.99
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
      var f = {
        callback: function() {}
      };
      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
      Ext.util.Observable.capture(ctrl.application, f.callback);

      var scans = [{
        "id": 1152,
        "rt": 1652.25
      }, {
        "id": 1133,
        "rt": 1624.99
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
          down: function() {
            return mocked_chromatogram;
          }
        };
      }
    };
    ctrl.center(tool);

    expect(mocked_chromatogram.resetScales).toHaveBeenCalled();
  });

  it('showUploadForm', function() {
    ctrl.hasStructures = false;
    spyOn(mocked_form_panel, 'setDisabledAnnotateFieldset');
    spyOn(mocked_form_panel, 'loadDefaults');
    var panel = {
      setActiveItem: function() {}
    };
    spyOn(panel, 'setActiveItem');
    spyOn(ctrl, 'getChromatogramPanel').andReturn(panel);

    ctrl.showUploadForm();

    expect(mocked_form_panel.setDisabledAnnotateFieldset).toHaveBeenCalledWith(true);
    expect(mocked_form_panel.loadDefaults).toHaveBeenCalledWith(appRootBase + '/data/runinfo.json');
    expect(panel.setActiveItem).toHaveBeenCalledWith(1);
  });

  it('showChromatogram', function() {
    var panel = {
      setActiveItem: function() {}
    };
    spyOn(panel, 'setActiveItem');
    spyOn(ctrl, 'getChromatogramPanel').andReturn(panel);

    ctrl.showChromatogram();

    expect(panel.setActiveItem).toHaveBeenCalledWith(0);
  });

  it('uploadHandler', function() {
    spyOn(mocked_form, 'isValid').andReturn(true);
    spyOn(mocked_form, 'submit');

    ctrl.uploadHandler();

    expect(mocked_form.isValid).toHaveBeenCalledWith();
    expect(mocked_form.submit).toHaveBeenCalledWith({
      url: appRootBase + '/rpc/' + ctrl.application.jobid + '/add_ms_data',
      submitEmptyText: false,
      waitMsg: jasmine.any(String),
      success: jasmine.any(Function),
      failure: jasmine.any(Function)
    });
  });

  it('uploadHandler with invalid form', function() {
    spyOn(mocked_form, 'isValid').andReturn(false);
    spyOn(mocked_form, 'submit');

    ctrl.uploadHandler();

    expect(mocked_form.isValid).toHaveBeenCalledWith();
    expect(mocked_form.submit).not.toHaveBeenCalled();
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

  describe('onZoomDirectionChange', function() {
    it('x', function() {
      spyOn(ctrl, 'onZoomDirectionChange');

      ctrl.onZoomDirectionXChange(null, false);

      expect(ctrl.onZoomDirectionChange).toHaveBeenCalledWith('x', false);
    });

    it('y', function() {
      spyOn(ctrl, 'onZoomDirectionChange');

      ctrl.onZoomDirectionYChange(null, false);

      expect(ctrl.onZoomDirectionChange).toHaveBeenCalledWith('y', false);
    });

    it('setZoom', function() {
      spyOn(mocked_chromatogram, 'setZoom');

      ctrl.onZoomDirectionChange('z', true);

      expect(mocked_chromatogram.setZoom).toHaveBeenCalledWith('z', true);
    });
  });

  describe('applyRole', function() {
    it('canrun', function() {
      ctrl.application.features.run = true;
      spyOn(ctrl.actionsMenu, 'hideUploadAction');

      ctrl.applyRole();

      expect(ctrl.actionsMenu.hideUploadAction).not.toHaveBeenCalledWith();
    });

    it('cantrun', function() {
      ctrl.application.features.run = false;
      spyOn(ctrl.actionsMenu, 'hideUploadAction');

      ctrl.applyRole();

      expect(ctrl.actionsMenu.hideUploadAction).toHaveBeenCalledWith();
    });
  });

  it('loadExample', function() {
    spyOn(ctrl.application, 'runInfoUrl').andReturn('http://example.com/defaults.json');
    spyOn(mocked_form, 'load');

    ctrl.loadExample();

    var expected = {
      url: 'http://example.com/defaults.json?selection=example',
      method: 'GET',
      waitMsg: 'Fetching example settings',
      failure: jasmine.any(Function)
    };
    expect(mocked_form.load).toHaveBeenCalledWith(expected);
  });

  describe('changeMsDataFormat', function() {
    var ms_data_format_field = null;

    beforeEach(function() {
      var names = [
        'filter_heading',
        'max_ms_level',
        'abs_peak_cutoff',
        'scan',
        'precision_heading',
        'mz_precision',
        'mz_precision_abs',
        'precursor_mz_precision',
        'intensity_heading',
        'ms_intensity_cutoff',
        'msms_intensity_cutoff'
      ];
      mocked_form.fields = new Ext.util.MixedCollection();
      var i = 0;
      Ext.Array.forEach(names, function(name) {
        var f = {
          name: name,
          getName: function() {
            return this.name;
          },
          disabled: false,
          hidden: false,
          enable: function() {
            this.disabled = false;
          },
          disable: function() {
            this.disabled = true;
          },
          show: function() {
            this.hidden = false;
          },
          hide: function() {
            this.hidden = true;
          }
        };
        mocked_form.fields.add(i++, f);
      });
      mocked_form.getFields = function() {
        return this.fields;
      };
      mocked_form.getFieldStates = function() {
        var state = {};
        this.fields.each(function(f) {
          state[f.getName()] = {
            hidden: f.hidden,
            disabled: f.disabled
          };
        });
        return state;
      };

      ms_data_format_field = {};
      ms_data_format_field.up = function() {
        return mocked_form_panel;
      };
    });

    it('form_tree', function() {
      ctrl.changeMsDataFormat(ms_data_format_field, 'form_tree');

      var estates = {
        filter_heading: {
          hidden: true,
          disabled: true
        },
        max_ms_level: {
          hidden: true,
          disabled: true
        },
        abs_peak_cutoff: {
          hidden: true,
          disabled: true
        },
        scan: {
          hidden: true,
          disabled: true
        },
        precision_heading: {
          hidden: true,
          disabled: true
        },
        mz_precision: {
          hidden: true,
          disabled: true
        },
        mz_precision_abs: {
          hidden: true,
          disabled: true
        },
        precursor_mz_precision: {
          hidden: true,
          disabled: true
        },
        intensity_heading: {
          hidden: true,
          disabled: true
        },
        ms_intensity_cutoff: {
          hidden: true,
          disabled: true
        },
        msms_intensity_cutoff: {
          hidden: true,
          disabled: true
        }
      };
      expect(mocked_form.getFieldStates()).toEqual(estates);
    });

    it('mass_tree', function() {
      ctrl.changeMsDataFormat(ms_data_format_field, 'mass_tree');

      var estates = {
        filter_heading: {
          hidden: true,
          disabled: true
        },
        max_ms_level: {
          hidden: true,
          disabled: true
        },
        abs_peak_cutoff: {
          hidden: true,
          disabled: true
        },
        scan: {
          hidden: true,
          disabled: true
        },
        precision_heading: {
          hidden: false,
          disabled: false
        },
        mz_precision: {
          hidden: false,
          disabled: false
        },
        mz_precision_abs: {
          hidden: false,
          disabled: false
        },
        precursor_mz_precision: {
          hidden: true,
          disabled: true
        },
        intensity_heading: {
          hidden: true,
          disabled: true
        },
        ms_intensity_cutoff: {
          hidden: true,
          disabled: true
        },
        msms_intensity_cutoff: {
          hidden: true,
          disabled: true
        }
      };
      expect(mocked_form.getFieldStates()).toEqual(estates);
    });

    it('mzxml', function() {
      ctrl.changeMsDataFormat(ms_data_format_field, 'form_tree');

      ctrl.changeMsDataFormat(ms_data_format_field, 'mzxml');

      var estates = {
        filter_heading: {
          hidden: false,
          disabled: false
        },
        max_ms_level: {
          hidden: false,
          disabled: false
        },
        abs_peak_cutoff: {
          hidden: false,
          disabled: false
        },
        scan: {
          hidden: false,
          disabled: false
        },
        precision_heading: {
          hidden: false,
          disabled: false
        },
        mz_precision: {
          hidden: false,
          disabled: false
        },
        mz_precision_abs: {
          hidden: false,
          disabled: false
        },
        precursor_mz_precision: {
          hidden: false,
          disabled: false
        },
        intensity_heading: {
          hidden: false,
          disabled: false
        },
        ms_intensity_cutoff: {
          hidden: false,
          disabled: false
        },
        msms_intensity_cutoff: {
          hidden: false,
          disabled: false
        }
      };
      expect(mocked_form.getFieldStates()).toEqual(estates);
    });
  });

  it('showHelp', function() {
    spyOn(ctrl.application, 'showHelp');

    ctrl.showHelp();

    expect(ctrl.application.showHelp).toHaveBeenCalledWith('chromatogram');
  });
});
