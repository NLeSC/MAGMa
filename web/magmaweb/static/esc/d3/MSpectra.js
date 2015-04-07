/**
 * A Mass spectra viewer
 * @class Esc.d3.MSpectra
 * @extends Esc.d3.Abstract
 * @author Stefan Verhoeven
 *
 * A example with 2 peaks, a marker at the first peak and select of first peak:
 *
 *     @example
 *     var mass_spectra = Ext.create('Esc.d3.MSpectra', {
 *       renderTo: Ext.getBody(),
 *       title: 'Mass spectra',
 *       width: 400, height: 300,
 *       axesPadding: [16, 5, 58, 80],
 *       data: [{mz:1,intensity:2},{mz:3,intensity:4}],
 *       markers: [{mz: 3}],
 *       cutoff: 1
 *     });
 *     mass_spectra.selectPeak(1);
 *
 * Note! This example requires d3.v3.js to be sourced.
 *
 */
Ext.define('Esc.d3.MSpectra', {
  extend: 'Esc.d3.Abstract',
  alias: 'widget.mspectra',
  initComponent: function() {
    var defConfig = {
        /**
         * @cfg {Array} data array of objects with mz and intensity and hashit properties.
         */
        /**
         * @property {Number} selectedPeak mz of selectedPeak.
         * When no peaks are selected then it is set to false.
         * @readonly
         */
        selectedpeak: false,
        /**
         * @cfg {Number} cutoff intensity under which peaks where disregarded
         */
        cutoff: 0,
        /**
         * @cfg {String} cutoffCls The CSS class applied to cutoff line.
         */
        cutoffCls: 'cutoffline',
        /**
         * @cfg {String} selectedPeakCls The CSS class applied to markers of a selected peak.
         */
        selectedPeakCls: 'selected',
        markers: [],
        chartWidth: 0,
        chartHeight: 0
        // TODO add cutoff and basepeakintensity so cutoff line can be drawn
    };

    Ext.applyIf(this, defConfig);

    this.callParent(arguments);

    this.addEvents(
      /**
       * @event selectpeak
       * Fires when user clicks on peak marker (triangle)
       * @param {Number} mz
       */
      'selectpeak',
      /**
       * @event unselectpeak
       * Fires when user clicks on a selected peak marker (triangle) which unselects it.
       * @param {Number} mz
       */
      'unselectpeak',
      /**
       * @event
       * Fires when mouse is moved over a vertical line of the scan
       * @param {Object} peak
       * @param {Number} peak.mz M/z of peak
       * @param {Number} peak.intensity Intensity of peak.
       */
      'mouseoverpeak'
    );
  },
  /**
   * @inheritdoc Esc.d3.Abstract#onZoom
   */
  onZoom: function() {
    this.callParent(arguments);
    var me = this;

    if (this.markers.length) {
      this.svg.selectAll("path.lowermarker")
        .attr("transform", function(d) { return "translate(" + me.scales.x(d.mz) + "," + (me.scales.y(0)+4) + ")"; });
      this.svg.selectAll("path.uppermarker")
        .attr("transform", function(d) { return "translate(" + me.scales.x(d.mz) + "," + (me.scales.y(me.ranges.y.max)-4) + ")"; });
    }

    this.svg.selectAll("line.mspeak")
        .attr("x1", function(d) { return me.scales.x(d.mz); })
        .attr("y1", function(d) { return me.scales.y(0); })
        .attr("x2", function(d) { return me.scales.x(d.mz); })
        .attr("y2", function(d) { return me.scales.y(d.intensity); })
    ;

    this.svg.selectAll("line."+this.cutoffCls)
		.attr('x1',0)
		.attr('x2',this.chartWidth)
		.attr('y1',this.scales.y(this.cutoff))
		.attr('y2',this.scales.y(this.cutoff))
    ;
  },
  initScales: function() {
    this.callParent(arguments);
    this.ranges.x.min = 0;
    this.ranges.x.max = d3.max(this.data, function(r) { return r.mz; });
    this.ranges.y.min = 0;
    this.ranges.y.max = d3.max(this.data, function(r) { return r.intensity; });
    this.scales.x = d3.scale.linear().domain([this.ranges.x.min, this.ranges.x.max]).range([0, this.chartWidth]);
    this.scales.y = d3.scale.linear().domain([this.ranges.y.min, this.ranges.y.max]).range([this.chartHeight, 0]);
  },
  initAxes: function() {
	this.callParent(arguments);
    var nrxticks = this.ticks.x;
    if (this.chartWidth < 25*(6+2)) {
        nrxticks = 2;
    }
    if (this.chartWidth < 25*2) {
        nrxticks = 0;
    }
    this.axes.x = d3.svg.axis().scale(this.scales.x).ticks(nrxticks);

    var nryticks = this.ticks.y;
    if (this.chartHeight < 16*(1+2)) {
        nryticks = 2;
    }
    if (this.chartHeight < 16*2) {
        nryticks = 0;
    }
    this.axes.y = d3.svg.axis().scale(this.scales.y).ticks(nryticks).orient("left").tickFormat(d3.format('.2e'));
  },
  draw: function() {
    this.callParent(arguments);
    var me = this;

    // Add the x-axis.
    this.svg.append("svg:g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + this.chartHeight + ")")
        .call(this.axes.x)
        .append("svg:text")
          .attr("x",this.chartWidth/2).attr("y",30)
          .attr("text-anchor","middle")
          .text('M/z')
    ;

    // Add the y-axis.
    this.svg.append("svg:g")
        .attr("class", "y axis")
        .call(this.axes.y)
    ;

    // cutoff
    if (this.cutoff) {
      this.svg.append("svg:line")
        .attr('class', this.cutoffCls)
        .attr('x1',0)
        .attr('x2',this.chartWidth)
        .attr('y1',this.scales.y(this.cutoff))
        .attr('y2',this.scales.y(this.cutoff))
      ;
    }

    // of each mz plot intensity as vertical line
    var peaks = this.svg.selectAll('line.mspeak').data(this.data, function(d) { return d.mz;});

    peaks.enter().append("svg:line")
        .attr("class", "mspeak")
        .classed('assigned', function(d) { return d.assigned_molid !== null;})
        .attr("x1", function(d) { return me.scales.x(d.mz); })
        .attr("y1", function(d) { return me.scales.y(0); })
        .attr("x2", function(d) { return me.scales.x(d.mz); })
        .attr("y2", function(d) { return me.scales.y(d.intensity); })
        .on('mouseover', function(peak) {
            me.fireEvent('mouseoverpeak', peak);
        })
    ;

    if (this.hasMarkers()) {
      this.onMarkersReady();
    }
  },
  undraw: function() {
    this.svg.selectAll('g.axis').remove();
    this.svg.selectAll('line.mspeak').remove();
    this.clearPeakSelection();
    this.svg.selectAll('.marker').remove();
    this.svg.selectAll('.'+this.cutoffCls).remove();
    this.callParent(arguments);
  },
  /**
   * @private
   */
  onToggleMarker: function(mz) {
    var me = this;
    this.markerSelect(function(e) {
      return (mz == e.mz && me.selectedpeak != e.mz);
    });
    if (mz != me.selectedpeak) {
      if (me.selectedpeak) {
    	me.fireEvent('unselectpeak', me.selectedpeak);
      }
      me.fireEvent('selectpeak', mz);
      me.selectedpeak = mz;
    } else {
      me.fireEvent('unselectpeak', mz);
      me.selectedpeak = false;
    }
  },
  /**
   * selects markers based on f returning true or false.
   * @param f function
   * @private
   */
  markerSelect: function(f) {
    this.svg.selectAll("path.marker")
      .classed(this.selectedPeakCls, f)
    ;
  },
  /**
   * Select peak by its mz
   * @param mz float An mz of a peak
   */
  selectPeak: function(mz) {
	var sameMz = function(d) {
      return (d.mz == mz);
    };
    if (this.data.some(sameMz)) {
    	// only select peak when mspectra has peak with that mz
        this.markerSelect(sameMz);
    	this.selectedpeak = mz;
    	return true;
    } else {
    	return false;
    }
  },
  /**
   * Clears any selected peaks
   */
  clearPeakSelection: function() {
	var selectedPeak = this.selectedpeak;
    this.markerSelect(false);
    this.selectedpeak = false;
    if (selectedPeak) {
  	  this.fireEvent('unselectpeak', selectedPeak);
  	}
  },
  /**
   * Set markers on mz's which can be selected
   * @param data array of objects with mz prop
   */
  setMarkers: function(data) {
    this.clearPeakSelection();
    this.svg.selectAll('.marker').remove();
    this.markers = data;
    this.onMarkersReady();
  },
  hasMarkers: function() {
    return (this.markers.length>0);
  },
  onMarkersReady: function() {
    // can not add markers if there is no data
    if (!this.hasData()) {
      return;
    }

    var me = this;
    function markerTitle(d) {
      return 'm/z='+d.mz;
    }
    function markerClick(d) {
      me.onToggleMarker(d.mz);
    }

    // lower markers
    this.svg.selectAll("path.lowermarker")
    .data(function() {return me.markers;})
    .enter().append("svg:path")
      .attr('class', 'marker lowermarker annotated')
      .attr("transform", function(d) { return "translate(" + me.scales.x(d.mz) + "," + (me.scales.y(0)+4) + ")"; })
      .attr("d", d3.svg.symbol().type('triangle-up').size(32) )
      .on('click', markerClick)
      .append("svg:title")
        .text(markerTitle)
    ;

    // upper markers
    this.svg.selectAll("path.uppermarker")
    .data(function() {return me.markers;})
    .enter().append("svg:path")
      .attr('class', 'marker uppermarker annotated')
      .attr("transform", function(d) { return "translate(" + me.scales.x(d.mz) + "," + (me.scales.y(me.ranges.y.max)-4) + ")"; })
      .attr("d", d3.svg.symbol().type('triangle-down').size(32) )
      .on('click', markerClick)
      .append("svg:title")
        .text(markerTitle)
    ;
  }
});