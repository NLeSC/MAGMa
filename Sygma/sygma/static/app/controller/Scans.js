Ext.define('Esc.msygma.controller.Scans', {
  extend: 'Ext.app.Controller',
  views: [ 'scan.Chromatogram' ],
  refs: [{
    ref: 'scanChromatogram', selector: 'scanchromatogram'
  }],
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
  },
  onLaunch: function() {
    var me = this;
    console.log('Scans contr launch');
    // config chromatogram,
    // must be done after viewport so canvas is avaliable
    var chromatogram = this.getScanChromatogram();
    chromatogram.cutoff = this.application.getMs_intensity_cutoff();
    chromatogram.setLoading(true);
    d3.json(me.application.getUrls().chromatogram, function(data) {
      chromatogram.setLoading(false);
      console.log('Loading chromatogram');
      chromatogram.setData(data);
      me.resetScans();
    });
  },
  clearExtractedIonChromatogram: function() {
    this.getScanChromatogram().setExtractedIonChromatogram([]);
  },
  loadExtractedIonChromatogram: function(metid, metabolite) {
    console.log('Loading extracted ion chromatogram');
    this.getScanChromatogram().setLoading(true);
    var me = this;
    Ext.Ajax.request({
      url: Ext.String.format(this.application.getUrls().extractedionchromatogram, metid),
      success: function(response) {
        me.getScanChromatogram().setLoading(false);
        var obj = Ext.decode(response.responseText);
        me.getScanChromatogram().setExtractedIonChromatogram(obj.chromatogram);
        me.setScans(obj.scans);
      }
    });
  },
  searchScan: function() {
    var me = this;
    Ext.MessageBox.prompt('Scan#', 'Please enter a level 1 scan identifier:', function(b,v) {
      if (b != 'cancel' && v) {
       v = v*1;
       me.selectScan(v);
      }
    });
  },
  clearScanSelection: function() {
    this.getScanChromatogram().clearScanSelection();
    this.onUnselectScan();
  },
  selectScan: function(scanid) {
    this.getScanChromatogram().selectScans([scanid]);
    this.application.fireEvent('selectscan', scanid);
  },
  /**
   * Each time the metabolite grid is loaded the response also contains a list of scans
   * where the filtered metabolites have hits, we use this to mark the scans that can be selected
   */
  setScansOfMetabolites: function(metabolitestore) {
      this.scans_of_metabolites = metabolitestore.getProxy().getReader().rawData.scans;
      this.setScans(this.scans_of_metabolites);
  },
  /**
   * Sets scans markers to scans where current metabolite filter has hits.
   */
  resetScans: function() {
      if (this.scans_of_metabolites) {
          this.setScans(this.scans_of_metabolites);
      }
  },
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
         chromatogram.selectScans([selectedScan]);
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
