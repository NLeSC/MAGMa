/**
 * Scans controller.
 *
 * Handles actions performed on the scan views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.controller.Scans', {
  extend: 'Ext.app.Controller',
  views: [ 'scan.Panel' ],
  refs: [{
      ref: 'chromatogramPanel', selector: 'scanpanel'
  }, {
      ref: 'chromatogram', selector: 'chromatogram'
  }, {
      ref: 'uploadForm', selector: 'scanuploadform'
  }],
  uses: [
         'Ext.window.MessageBox',
         'Ext.window.Window',
         'Ext.form.Panel',
         'Esc.magmaweb.view.scan.UploadFieldSet',
         'Esc.magmaweb.view.fragment.AnnotateFieldSet'
  ],
  /**
   * Cached scans which belong to (filtered) metabolites
   */
  scans_of_metabolites: [],
  init: function() {
    Ext.log({}, 'Scans controller init');
    var me = this;

    this.control({
      'chromatogram': {
        selectscan: function(scanid) {
            me.application.fireEvent('selectscan', scanid);
        },
        unselectscan: function(scanid) {
            me.application.fireEvent('noselectscan', scanid);
        },
        mouseoverscan: function(scan) {
            if ('metaboliteintensity' in scan) {
              this.getChromatogramPanel().setTitle(Ext.String.format(
                      'Chromatogram (rt={0}, basepeak intensity={1}, metabolite intensity={2}, scan={3})',
                      scan.rt, scan.intensity, scan.metaboliteintensity, scan.id
              ));
            } else {
              this.getChromatogramPanel().setTitle(Ext.String.format('Chromatogram (rt={0}, intensity={1}, scan={2})', scan.rt, scan.intensity, scan.id));
            }
        }
      },
      'scanpanel tool[action=search]': {
        click: this.searchScan
      },
      'scanpanel tool[action=clearselection]': {
        click: this.clearScanSelection
      },
      'scanpanel tool[action=center]': {
        click: this.center
      },
      'scanpanel tool[action=actions]': {
        click: this.showActionsMenu
      },
      'scanuploadform component[action=uploadmsdata]': {
        click: this.uploadHandler
      },
      'scanuploadform component[action=uploadmsdatacancel]': {
        click: this.showChromatogram
      },
      'scanuploadform component[action=loadmsdataexample]': {
        click: this.loadExample
      },
	  'scanuploadform component[name=ms_data_format]': {
	    change: this.changeMsDataFormat
	  }
    });

    this.application.on('metaboliteload', this.setScansOfMetabolites, this);
    this.application.on('metaboliteselect', this.loadExtractedIonChromatogram, this);
    this.application.on('metabolitedeselect', function() {
        this.resetScans();
        this.clearExtractedIonChromatogram();
    }, this);
    this.application.on('metabolitenoselect', function() {
        this.resetScans();
        this.clearExtractedIonChromatogram();
    }, this);
    this.application.on('rpcsubmitsuccess', function() {
        Ext.getCmp('uploadmssaction').disable();
    });
    this.application.on('assignmentchanged', function(isAssigned, params) {
        me.loadChromatogram(function(data) {
            var chromatogram = this.getChromatogram();
            var selectedScan = chromatogram.selectedScan;
            chromatogram.setLoading(false);
            chromatogram.setData(data.scans);
            chromatogram.selectScan(selectedScan);
        });
    });

    this.actionsMenu = Ext.create('Ext.menu.Menu', {
        items: [{
            iconCls: 'icon-add',
            id: 'uploadmssaction',
            text: 'Upload MS data',
            handler: this.showUploadForm.bind(this)
        }, {
            text: 'Zoom direction',
            menu: {
            	items: [{
            		text: 'X axis',
            		checked: true,
            		checkHandler: this.onZoomDirectionXChange.bind(this)
            	}, {
            		text: 'Y axis',
            		checked: false,
            		checkHandler: this.onZoomDirectionYChange.bind(this)
            	}]
            }
        }],
        hideUploadAction: function() {
        	this.getComponent('uploadmssaction').hide();
        }
    });

    /**
     * @property {Boolean} hasStructures
     * Whether there are structures.
     * Used to disable/enable annotate options
     */
    this.hasStructures = false;

    this.application.addEvents(
      /**
       * @event
       * Triggered when a scan is selected.
       * @param {Number} scanid Scan identifier.
       */
      'selectscan',
      /**
       * @event
       * Triggered when no scan is selected anymore.
       */
      'noselectscan',
      /**
       * @event chromatogramload
       * Fired after chromatogram has been loaded
       * @param chromatogram Chromatogram
       */
      'chromatogramload'
    );

  },
  /**
   * Loads the chromatogram from server.
   */
  onLaunch: function() {
      this.loadChromatogram(this.loadChromatogramCallback);
      this.applyRole();
  },
  loadChromatogram: function(callback) {
      var me = this;
      // config chromatogram,
      // must be done after viewport so canvas is avaliable
      var chromatogram = this.getChromatogram();
      chromatogram.setLoading(true);
      d3.json(
          me.application.getUrls().chromatogram,
          callback.bind(me)
      );
  },
  /**
   * Callback for loading chromatogram
   *
   * @params {Array} data Array of scans with rt, id and itensity props
   */
  loadChromatogramCallback: function(data) {
    if (data === null) {
      Ext.Error.raise({
         msg: 'Failed to load chromatogram from server'
      });
      return false;
    }
    var me = this;
    var chromatogram = this.getChromatogram();
    chromatogram.setLoading(false);
    Ext.log({}, 'Loading chromatogram');
    if (data.cutoff !== null) {
        chromatogram.cutoff = data.cutoff;
    }
    chromatogram.setData(data.scans);
    me.resetScans();
    if (data.scans.length === 0) {
    	// when there are no scans then user should upload some
        this.showUploadForm();
    } else if (data.scans.length === 1) {
    	// hide when chromatogram consists of 1 scan
    	// on molecule load the scan will be selected
    	this.getChromatogramPanel().hide();
    }
    this.application.fireEvent('chromatogramload', chromatogram);
  },
  clearExtractedIonChromatogram: function() {
    this.getChromatogram().setExtractedIonChromatogram([]);
  },
  /**
   * Download the extracted ion chromatogram of a metabolite on the chromatogram.
   * @param {Number} metid Metobolite identifier
   */
  loadExtractedIonChromatogram: function(metid) {
    Ext.log({}, 'Loading extracted ion chromatogram');
    this.getChromatogram().setLoading(true);
    var me = this;
    d3.json(
      Ext.String.format(this.application.getUrls().extractedionchromatogram, metid),
      me.loadExtractedIonChromatogramCallback.bind(me)
    );
  },
  /**
   * Callback for loading a extracted ion chromatogram
   * @param {Object} data
   * @param {Array} data.scans Scans in which metabolite has hits
   * @param {Array} data.chromatogram Foreach scan the intensity of the peak with metabolite m/z
   */
  loadExtractedIonChromatogramCallback: function(data) {
    if (data === null) {
      Ext.Error.raise({
         msg: 'Failed to load extracted ion chromatogram from server'
      });
      return false;
    }
    var chromatogram = this.getChromatogram();
    chromatogram.setLoading(false);
    chromatogram.setExtractedIonChromatogram(data.chromatogram);
    this.setScans(data.scans);
  },
  searchScan: function() {
    var me = this;
    Ext.MessageBox.prompt(
      'Scan#',
      'Please enter a level 1 scan identifier:',
      function(b,v) {
        if (b != 'cancel' && v) {
         v = v*1;
         me.selectScan(v, true);
        }
      }
    );
  },
  clearScanSelection: function() {
    this.getChromatogram().clearScanSelection();
    this.application.fireEvent('noselectscan');
  },
  /**
   * Select a scan
   * @param {Number} scanid Scan identifier
   * @param {Boolean} [silent=false] Passing true will supress the 'unselectscan' event from being fired.
   */
  selectScan: function(scanid, silent) {
    this.getChromatogram().selectScan(scanid, silent);
    this.application.fireEvent('selectscan', scanid);
  },
  /**
   * Each time the metabolite grid is loaded the response also contains a list of scans
   * where the filtered metabolites have hits, we use this to mark the scans that can be selected
   *
   * @param {Esc.magmaweb.store.Metabolites} metabolitestore rawdata of store reader has scans
   */
  setScansOfMetabolites: function(metabolitestore) {
      this.hasStructures = metabolitestore.getTotalCount() > 0;
      this.scans_of_metabolites = metabolitestore.getProxy().getReader().rawData.scans;
      this.setScans(this.scans_of_metabolites);
  },
  /**
   * Sets scans markers to scans where current metabolite filter has hits.
   */
  resetScans: function() {
    this.setScans(this.scans_of_metabolites);
  },
  /**
   * Add scan markers to chromatogram that can be selected.
   * @param {Array} scans Array of scans
   */
  setScans: function(scans) {
    Ext.log({}, 'Setting chromatogram scan markers');
    var chromatogram = this.getChromatogram();
    if (!chromatogram.hasData()) {
        return; // can not set scan markers if chromatogram is not loaded
    }
    if (scans.length) {
       // if scan is already selected and is part of new scans then reselect scan
       if (
         scans.some(function(e) {
           return (e.id == chromatogram.selectedScan);
         })
       ) {
         var selectedScan = chromatogram.selectedScan;
         chromatogram.setMarkers(scans);
         chromatogram.selectScan(selectedScan);
       } else {
         chromatogram.setMarkers(scans);
       }
       // if only one scan then select scan
       if (scans.length == 1 && chromatogram.selectedScan != scans[0].id) {
         this.selectScan(scans[0].id);
       }
    } else {
      this.application.fireEvent('noscansfound');
    }
  },
  center: function() {
      this.getChromatogram().resetScales();
  },
  showUploadForm: function() {
      var me = this;
      this.getUploadForm().setDisabledAnnotateFieldset(!this.hasStructures);
      this.getUploadForm().loadDefaults(me.application.runInfoUrl());
      this.getChromatogramPanel().setActiveItem(1);
  },
  showChromatogram: function() {
      this.getChromatogramPanel().setActiveItem(0);
  },
  uploadHandler: function() {
      var me = this;
      var form = this.getUploadForm().getForm();
      if (form.isValid()) {
          form.submit({
              url: this.application.rpcUrl('add_ms_data'),
              waitMsg: 'Submitting action ...',
              submitEmptyText: false,
              success: function(fp, o) {
                  var response = Ext.JSON.decode(o.response.responseText);
                  me.application.fireEvent('rpcsubmitsuccess', response.jobid);
              },
              failure: function(form, action) {
                  if (action.failureType === "server") {
                      Ext.Error.raise(Ext.JSON.decode(action.response.responseText));
                    } else {
                      Ext.Error.raise(action.response.responseText);
                    }
              }
          });
      }
  },
  /**
   * Show actions menu at event xy
   * @param {Ext.Element} tool
   * @param {Ext.EventObject} event
   */
  showActionsMenu: function(tool, event) {
     this.actionsMenu.showAt(event.getXY());
  },
  /**
   * Enable/Disable zoom on X axis
   * @param {Ext.Component} item
   * @param {Boolean} checked
   */
  onZoomDirectionXChange: function(item, checked) {
	  this.onZoomDirectionChange('x', checked);
  },
  /**
   * Enable/Disable zoom on Y axis
   * @param {Ext.Component} item
   * @param {Boolean} checked
   */
  onZoomDirectionYChange: function(item, checked) {
	  this.onZoomDirectionChange('y', checked);
  },
  /**
   * Enable/Disable zoom on a axis
   * @param {String} axes Name of axis. Can be 'x' or 'y'.
   * @param {Boolean} checked
   */
  onZoomDirectionChange: function(axis, checked) {
	  this.getChromatogram().setZoom(axis, checked);
  },
  /**
   * Apply role to user interface.
   * Checks canRun and if false removes all action buttons.
   */
  applyRole: function() {
	  if (this.application.canRun) {
		  return;
	  }
	  this.actionsMenu.hideUploadAction();
	  // TODO change tooltip of gears tool
  },
  /**
   * In MS Data upload forms loads the example data set.
   */
  loadExample: function() {
	  var form = this.getUploadForm().getForm();
	  var example_url = this.application.runInfoUrl()+'?selection=example';
	  form.load({
		  url: example_url,
	      method: 'GET',
	      waitMsg: 'Fetching example settings',
	      failure: function(form, action) {
              Ext.Error.raise(action.response.responseText);
          }
	  });
  },
  /**
   * Called when MS data format is changed.
   * When the 'tree' format is chosen the filtering is disabled.
   * When a non-'tree' format is chosen filtering is left intact.
   */
  changeMsDataFormat: function(field, value) {
      var form = this.getUploadForm().getForm();
	  if (value == 'tree') {
		  form.setValues({
		      'ms_intensity_cutoff': 0,
		      'msms_intensity_cutoff': 0,
		      'abs_peak_cutoff': 0
		  });
	  }
  }
});
