/**
 * Scans controller.
 *
 * Handles actions performed on the scan views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.mmm.controller.Scans', {
  extend: 'Ext.app.Controller',
  views: [ 'scan.Chromatogram' ],
  refs: [{
    ref: 'scanChromatogram', selector: 'scanchromatogram'
  }],
  uses: [ 'Ext.window.MessageBox' ],
  /**
   * Cached scans which belong to (filtered) metabolites
   */
  scans_of_metabolites: [],
  init: function() {
    console.log('Scans controller init');
    var me = this;

    this.control({
      'scanchromatogram': {
        selectscan: function(scanid) {
            me.application.fireEvent('selectscan', scanid);
        },
        unselectscan: function(scanid) {
            me.application.fireEvent('noselectscan', scanid);
        },
        mouseoverscan: function(scan) {
            if ('metaboliteintensity' in scan) {
              this.getScanChromatogram().setTitle(Ext.String.format(
                      'Chromatogram (rt={0}, basepeak intensity={1}, metabolite intensity={2}, scan={3})',
                      scan.rt, scan.intensity, scan.metaboliteintensity, scan.id
              ));
            } else {
              this.getScanChromatogram().setTitle(Ext.String.format('Chromatogram (rt={0}, intensity={1}, scan={2})', scan.rt, scan.intensity, scan.id));
            }
        }
      },
      'scanchromatogram tool[action=search]': {
        click: this.searchScan
      },
      'scanchromatogram tool[action=clearselection]': {
        click: this.clearScanSelection
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
      'noselectscan'
    );
  },
  /**
   * Loads the chromatogram from server.
   */
  onLaunch: function() {
    var me = this;
    console.log('Scans contr launch');
    // config chromatogram,
    // must be done after viewport so canvas is avaliable
    var chromatogram = this.getScanChromatogram();
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
    var me = this;
    var chromatogram = this.getScanChromatogram();
    chromatogram.setLoading(false);
    console.log('Loading chromatogram');
    chromatogram.setData(data);
    me.resetScans();
  },
  clearExtractedIonChromatogram: function() {
    this.getScanChromatogram().setExtractedIonChromatogram([]);
  },
  /**
   * Download the extracted ion chromatogram of a metabolite on the chromatogram.
   * @param {Number} metid Metobolite identifier
   */
  loadExtractedIonChromatogram: function(metid) {
    console.log('Loading extracted ion chromatogram');
    this.getScanChromatogram().setLoading(true);
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
    this.getScanChromatogram().setLoading(false);
    this.getScanChromatogram().setExtractedIonChromatogram(data.chromatogram);
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
    this.getScanChromatogram().clearScanSelection();
    this.application.fireEvent('noselectscan');
  },
  /**
   * Select a scan
   * @param {Number} scanid Scan identifier
   */
  selectScan: function(scanid) {
    this.getScanChromatogram().selectScan(scanid);
    this.application.fireEvent('selectscan', scanid);
  },
  /**
   * Each time the metabolite grid is loaded the response also contains a list of scans
   * where the filtered metabolites have hits, we use this to mark the scans that can be selected
   *
   * @param {Esc.mmm.store.Metabolites} metabolitestore rawdata of store reader has scans
   */
  setScansOfMetabolites: function(metabolitestore) {
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
    var chromatogram = this.getScanChromatogram();
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
  }
});
