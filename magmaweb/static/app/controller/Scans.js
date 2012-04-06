/**
 * Scans controller.
 *
 * Handles actions performed on the scan views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.controller.Scans', {
  extend: 'Ext.app.Controller',
  views: [ 'scan.Chromatogram' ],
  refs: [{
      ref: 'chromatogramPanel', selector: 'chromatogrampanel'
  }, {
      ref: 'chromatogram', selector: 'chromatogram'
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
    console.log('Scans controller init');
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
      'chromatogrampanel tool[action=search]': {
        click: this.searchScan
      },
      'chromatogrampanel tool[action=clearselection]': {
        click: this.clearScanSelection
      },
      'chromatogrampanel tool[action=center]': {
        click: this.center
      },
      'chromatogrampanel tool[action=upload]': {
        click: this.showUploadForm
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
    var me = this;
    // config chromatogram,
    // must be done after viewport so canvas is avaliable
    var chromatogram = this.getChromatogram();
    chromatogram.cutoff = this.application.getMs_intensity_cutoff();
    chromatogram.setLoading(true);
    d3.json(
        me.application.getUrls().chromatogram,
        me.loadChromatogramCallback.bind(me)
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
    console.log('Loading chromatogram');
    chromatogram.setData(data);
    me.resetScans();
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
    console.log('Loading extracted ion chromatogram');
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
         me.selectScan(v);
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
   */
  selectScan: function(scanid) {
    this.getChromatogram().selectScan(scanid);
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
    console.log('Setting chromatogram scan markers');
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
      if (!this.uploadForm) {
          this.uploadForm = Ext.create('Ext.window.Window', {
              title: 'Upload MS data',
              height: 320,
              width: 600,
              layout: 'fit',
              modal: true,
              closeAction: 'hide',
              items: {
                  xtype: 'form',
                  bodyPadding: 5,
                  defaults: { bodyPadding: 5 },
                  border: false,
                  autoScroll: true,
                  url: me.rpcUrl('add_ms_data'),
                  items: [{
                      xtype: 'uploadmsdatafieldset'
                  }, {
                      xtype : 'annotatefieldset',
                      disabled: !this.hasStructures,
                      collapsed : true,
                      collapsible : true
                  }],
                  buttons: [{
                      text: 'Submit',
                      handler: function() {
                          var form = this.up('form').getForm();
                          var wf = this.up('window');
                          if (form.isValid()) {
                              form.submit({
                                  waitMsg: 'Submitting action ...',
                                  success: function(fp, o) {
                                      console.log('Action submitted');
                                      wf.hide();
                                  },
                                  failure: function(form, action) {
                                      console.log(action.failureType);
                                      console.log(action.result);
                                      wf.hide();
                                  }
                              });
                          }
                      }
                  }, {
                      text: 'Reset',
                      handler: function() {
                          this.up('form').getForm().reset();
                      }
                  }]
              }
          });
      }
      this.uploadForm.show();
  }
});
