/**
 * MSpectras controller.
 *
 * Handles actions performed on the mspectra views.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.controller.MSpectras', {
  extend: 'Ext.app.Controller',
  requires: [ 'Esc.d3.MSpectra' ],
  config: {
    /**
     * Maximum MS level or nr of MS levels.
     * @cfg {Number}
     */
    maxmslevel: null,
    /**
     * MSpectra endpoint.
     * Tokenized string with scanid and mslevel tokens.
     * @cfg {String}
     */
    url: null,
  },
  /**
   * @property {Array} mspectras Array of Ext.esc.MSpectra
   * Index is MS level.
   */
  mspectras: [],
  constructor: function(config) {
    this.initConfig(config);
    this.callParent(arguments);
    return this;
  },
  init: function() {
    this.setUrl(this.application.getUrls().mspectra);
    this.setMaxmslevel(this.application.getMaxmslevel());
    this.initMSpectras();

    this.application.on('selectscan', this.loadMSpectra1, this);
    this.application.on('fragmentexpand', this.loadMSpectraFromFragment, this);
    this.application.on('fragmentload', this.loadMSpectrasFromFragment, this);
    this.application.on('fragmentselect', this.selectPeakFromFragment, this);
    this.application.on('fragmentdeselect', this.deselectPeakFromFragment, this);
    this.application.on('noselectscan', this.clearMSpectra1, this);
    this.application.on('scanandmetabolitenoselect', this.clearMSpectraFrom2, this);
    this.application.on('fragmentcollapse', this.clearMSpectraFromFragment, this);

    this.addEvents(
      /**
       * @event
       * Triggered when a peak in a MSpectra is selected.
       * @param {Number} mz M/z of selected peak
       * @param {Number} mslevel MS level where peak is located
       */
      'peakselect',
      /**
       * @event
       * Triggered when a peak in a MSpectra is deselected.
       * @param {Number} mz M/z of selected peak
       * @param {Number} mslevel MS level where peak is located
       */
      'peakdeselect',
      /**
       * @event
       * Triggered when a mspectra is loaded from the server
       * @param {Number} scanid Scan identifier
       * @param {Number} mslevel MS level of scan
       */
      'mspectraload',
      /**
       * @event
       * Triggered when a mspectra is cleared
       * @param {Number} mslevel MS level of scan
       */
      'mspectraclear',
      /**
       * @event
       * Fires when mouse is moved over a vertical line of a mspectra
       * @param {Object} peak
       * @param {Number} peak.mz M/z of peak
       * @param {Number} peak.intensity Intensity of peak.
       * @param {Number} mslevel MS level of scan
       * @param {Number} scanid Scan identifier of mspectra which is loaded
       */
      'peakmouseover'
    );

    // register controls foreach mspectra
    for (var mslevel = 1; mslevel <= this.getMaxmslevel(); mslevel++) {
        var centerquery = '#mspectra'+mslevel+'panel tool[action=center]';
        this.control(centerquery, { click: this.center });
    }
  },
  /**
   * Initializes MSpectra views
   * @private
   */
  initMSpectras: function() {
    var app = this.application;
    for (var mslevel = 1; mslevel <= this.getMaxmslevel(); mslevel++) {
      this.mspectras[mslevel] = Ext.create('Esc.d3.MSpectra', {
        mslevel: mslevel,
        emptyText: (
            mslevel==1 ?
            'Select a scan in the chromatogram' :
            'Select a fragment to show its level '+mslevel+' scan'
        ),
        listeners: {
          selectpeak: function(mz) {
            app.fireEvent('peakselect', mz, this.mslevel);
            // TODO unselect peaks of child scans
          },
          unselectpeak: function(mz) {
            app.fireEvent('peakdeselect', mz, this.mslevel);
          },
          mouseoverpeak: function(peak) {
            app.fireEvent('peakmouseover', peak, this.mslevel, this.scanid);
          }
        }
      });
    }
  },
  /**
   * Return MSpectra view based on MS level
   * @param {Number} mslevel
   */
  getMSpectra: function(mslevel) {
    return this.mspectras[mslevel];
  },
  /**
   * Load a MSpectra.
   *
   * @param {Number} mslevel MS level
   * @param {Number} scanid Scan identifier.
   * @param {Array} markers Array of markers to add after Mspectra is loaded.
   */
  loadMSpectra: function(mslevel, scanid, markers) {
    var me = this;
    console.log('Loading msspectra level '+mslevel+' with id '+scanid);
    this.getMSpectra(mslevel).setLoading(true);
    d3.json(
      Ext.String.format(this.getUrl(), scanid, mslevel),
      function(data) {
        me.onLoadMSpectra(mslevel, scanid, markers, data);
      }
    );
  },
  /**
   * onLoad a MSpectra.
   *
   * @param {Number} mslevel MS level
   * @param {Number} scanid Scan identifier.
   * @param {Array} markers Array of markers to add after Mspectra is loaded.
   * @param {Object} data
   * @param {Array} data.peaks Array of peaks of scan
   * @param {Number} data.cutoff Cutoff of this scan (based on basepeak intensity and msms intensity cutoff ratio)
   */
  onLoadMSpectra: function(mslevel, scanid, markers, data) {
    var mspectra = this.getMSpectra(mslevel);
    if (!data) {
      Ext.Error.raise({
          msg: 'Unable to find mspectra scan on level '+mslevel+' with id '+scanid,
      });
      return;
    }
    mspectra.setLoading(false);
    mspectra.scanid = scanid;
    mspectra.cutoff = data.cutoff;
    mspectra.setData(data.peaks);
    mspectra.setMarkers(markers);
    mspectra.up('panel').down('tool[action=center]').enable();
    this.application.fireEvent('mspectraload', scanid, mslevel);
  },
  /**
   * Load a lvl1 MSpectra
   *
   * @param {Number} scanid Scan identifier.
   */
  loadMSpectra1: function(scanid) {
    this.loadMSpectra(1, scanid, []);
  },
  /**
   * Load a lvl2 MSpectra
   *
   * @param {Number} scanid Scan identifier.
   * @param {Array} markers Array of markers to add after Mspectra is loaded.
   */
  loadMSpectra2: function(scanid, markers) {
    this.loadMSpectra(2, scanid, markers);
  },
  /**
   * Loads a MSpectra of a fragment.
   *
   * @param {Esc.magmaweb.model.Fragment} fragment MSpectra of fragments firstchild  is loaded
   */
  loadMSpectraFromFragment: function(fragment) {
    var mslevel = fragment.firstChild.data.mslevel;
    var scanid = fragment.firstChild.data.scanid;
    var mspectra = this.getMSpectra(mslevel);
    // skip if mspectra already has same scan
    if (mspectra.scanid != scanid) {
      this.loadMSpectra(
        mslevel, scanid,
        fragment.childNodes.map(function(r) { return {mz: r.data.mz}; })
      );
    }
  },
  /**
   * Depending on which mslevel the fragment was found and if has children
   * MSpectra are loaded and peaks selected.
   *
   * @param {Esc.magmaweb.model.Fragment} parent
   * @param {Array} children child fragments
   */
  loadMSpectrasFromFragment: function(parent, children) {
    // TODO remove ... || when extjs 4.1 is final
    if (('id' in parent.data && parent.data.id == 'root' ) || parent.data.root) {
      // lvl1 fragment
      console.log('Selecting metabolite peak in lvl1 mspectra and loading lvl2 mspectra');
      var mspectra = this.getMSpectra(1);
      // set markers in lvl1 scan
      mspectra.setMarkers(
        children.map(function(r) { return {mz: r.data.mz}; })
      );
      // lvl1 fragment is a peak in lvl1 scan
      var metabolite_fragment = parent.firstChild;
      mspectra.selectPeak(metabolite_fragment.data.mz);
      // with a optional child lvl2 scan
      if (metabolite_fragment.hasChildNodes()) {
        this.loadMSpectra2(
            metabolite_fragment.firstChild.data.scanid,
            metabolite_fragment.childNodes.map(
                function(r) { return {mz: r.data.mz}; }
            )
        );
      }
    } else {
      // lvl >1 fragment
      console.log('Select peak of fragment');
      this.getMSpectra(parent.data.mslevel).selectPeak(parent.data.mz);
    }
  },
  /**
   * Clear a MSpesctra.
   *
   * @param {Number} mslevel Level of MSpectra to clear
   */
  clearMSpectra: function(mslevel) {
    var mspectra = this.getMSpectra(mslevel);
    mspectra.setData([]);
    mspectra.scanid = -1;
    this.application.fireEvent('mspectraclear', mslevel);
    mspectra.up('panel').down('tool[action=center]').disable();
  },
  /**
   * Clear MSpectra lvl 1
   */
  clearMSpectra1: function() {
    this.clearMSpectra(1);
  },
  /**
   * @param {Number} mslevel This mspectra lvl is clear plus any higher lvl mspectra
   * @private
   */
  clearMSpectraFrom: function(mslevel) {
    for (var i = mslevel; i <= this.getMaxmslevel(); i++) {
      this.clearMSpectra(i);
    }
  },
  /**
   * Clears all mspectra >=lvl2 and clears lvl1 markers
   */
  clearMSpectraFrom2: function() {
    this.getMSpectra(1).setMarkers([]);
    this.clearMSpectraFrom(2);
  },
  /**
   * When fragment is collapsed the mspectra of its children must be cleared.
   * @param {Esc.magmaweb.model.Fragment} fragment
   */
  clearMSpectraFromFragment: function(fragment) {
    this.clearMSpectraFrom(fragment.data.mslevel+1);
  },
  /**
   * Select a peak using a fragment.
   * Peaks of fragments parents are also selected.
   * And peaks in child mspectras are unselected
   *
   * @param {Esc.magmaweb.model.Fragment} fragment Fragment of peak to select.
   */
  selectPeakFromFragment: function(fragment) {
    var mslevel = fragment.data.mslevel;
    this.getMSpectra(mslevel).selectPeak(fragment.data.mz);
    // clear selection in child mspectra
    for (var i = mslevel+1; i <= this.getMaxmslevel(); i++) {
      this.getMSpectra(i).clearPeakSelection();
    }
    // select parent peaks
    if (mslevel == 2 ) {
      this.getMSpectra(1).selectPeak(fragment.parentNode.data.mz);
    } else if (mslevel == 3) {
      this.getMSpectra(2).selectPeak(fragment.parentNode.data.mz);
      this.getMSpectra(1).selectPeak(fragment.parentNode.parentNode.data.mz);
    } else if (mslevel > 3) {
      // TODO use recursive func to select parent peaks
    }
  },
  /**
   * Deselect a peak using a fragment.
   *
   * @param {Esc.magmaweb.model.Fragment} fragment Fragment of peak to deselect.
   */
  deselectPeakFromFragment: function(fragment) {
    this.getMSpectra(fragment.data.mslevel).clearPeakSelection();
  },
  center: function(tool) {
      var mspectra = tool.up('panel').down('mspectra');
      mspectra.resetScales();
  }
});