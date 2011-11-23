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
        selectscan: this.onSelectScan,
        unselectscan: this.onUnselectScan
      },
      'scanchromatogram tool[action=search]': {
        click: this.searchScan
      },
      'scanchromatogram tool[action=clearselection]': {
        click: this.clearScanSelection
      }
    });
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
      me.setChromatogramMarkersByMetaboliteFilter();
    });
  },
  setChromatogramMarkersByMetaboliteFilter: function() {
    var store = this.application.getController('Metabolites').getMetabolitesStore();
    if (store.isLoaded && this.getScanChromatogram().hasData()) {
      console.log('Setting chromatogram markers');
      var markers = store.getProxy().getReader().rawData.scans;
      this.getScanChromatogram().setMarkers(markers);
    }
  },
  clearExtractedIonChromatogram: function() {
    this.getScanChromatogram().setExtractedIonChromatogram([]);
  },
  onSelectScan: function(scanid) {
    var me = this;
    console.log('select scan '+scanid);
    me.getController('Metabolites').getMetabolitesStore().removeScanFilter();
    me.getController('Fragments').clearFragments();
    // if metabolite has been selected
    // and scanid is hit of metabolite then show fragments with this metabolite and scan
    // else filter mstore and load scan
    if (
      me.getController('Metabolites').getMetaboliteList().getSelectionModel().selected.getCount() > 0
      &&
      me.getController('Metabolites').getMetaboliteList().getSelectionModel().selected.getAt(0).data.scans.some(
        function(e) { return (e.id == scanid); }
      )
    ) {
      me.application.loadMSpectra1(scanid, function() {
        me.getController('Fragments').loadFragments(
          scanid,
          me.getController('Metabolites').getMetaboliteList().getSelectionModel().selected.getAt(0).data.metid
        );
      });
    } else {
      me.application.loadMSpectra1(scanid, function() {
        me.getController('Metabolites').getMetabolitesStore().setScanFilter(scanid);
      });
      // TODO in spectra add markers for metabolites present in scan
    }
  },
  onUnselectScan: function() {
    var me = this;
    console.log('Unselect scan');
    me.getController('Metabolites').getMetabolitesStore().removeScanFilter();
    me.getController('Fragments').clearFragments();
    me.application.clearMSpectra(1);
  },
  searchScan: function() {
    var me = this;
    Ext.MessageBox.prompt('Scan#', 'Please enter a level 1 scan identifier:', function(b,v) {
      if (b != 'cancel' && v) {
       v = v*1;
       me.getScanChromatogram().selectScans([v]);
       me.onSelectScan(v);
      }
    });
  },
  clearScanSelection: function() {
     this.getScanChromatogram().clearScanSelection();
     this.onUnselectScan();
  }
});
