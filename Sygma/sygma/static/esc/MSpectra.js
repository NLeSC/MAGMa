/**
 * A Mass spectra viewer
 * @class Ext.esc.MSpectra
 * @extends Ext.Panel
 * @author Stefan Verhoeven
 */
Ext.define('Ext.esc.MSpectra', {
  extend: 'Ext.esc.AbstractD3',
  alias: 'widget.mspectra',
  initComponent: function() {
    var defConfig = {
        /**
         * @cfg {Array} data array of objects with mz and intensity and hashit properties.
         */
        // mz of selectedPeak
        selectedPeak: -1,
        /**
         * @cfg {String} selectedPeakCls The CSS class applied to markers of a selected peak.
         */
        /**
         * @cfg {Number} cutoff intensity under which peaks where disregarded
         */
        cutoff: 0,
        /**
         * @cfg {String} cutoffCls The CSS class applied to cutoff line.
         */
        cutoffCls: 'cutoffline',
        selectedPeakCls: 'selectedpeak',
        markers: [],
        chartWidth: 0,
        chartHeight: 0,
        // TODO add cutoff and basepeakintensity so cutoff line can be drawn
    };

    Ext.applyIf(this, defConfig);

    this.callParent(arguments);

    this.addEvents(
      /**
       * @event selectpeak
       * Fires when user clicks on peak marker (triangle)
       * @param {Int} mz
       */
      'selectpeak',
      /**
       * @event unselectpeak
       * Fires when user clicks on a selected peak marker (triangle) which unselects it.
       * @param {Int} mz
       */
      'unselectpeak'
    );
  },
  /**
   * @inheritdoc
   */
  redraw: function() {
    var me = this;
    if (d3.event && d3.event.translate[0] != 0 && d3.event.translate[1] != 0) {
      // pan and zoom x axis
      d3.event.transform(this.scales.x);
    }
    this.svg.select(".x.axis").call(this.axes.x);
    // do not scale y axis
    //svg.select(".y.axis").call(yAxis);
    if (this.markers.length) {
      this.svg.selectAll("path.lowermarker")
        .attr("transform", function(d) { return "translate(" + me.scales.x(d.mz) + "," + me.scales.y(0) + ")"; });
      this.svg.selectAll("path.uppermarker")
       .attr("transform", function(d) { return "translate(" + me.scales.x(d.mz) + "," + me.scales.y(me.ranges.y.max) + ")"; });
    }
    this.svg.selectAll("line.mspeak").attr("y2", function(d) {
      return me.scales.y(d.intensity);
    }).attr("x1", function(d) {
      return me.scales.x(d.mz);
    }).attr("x2", function(d) {
      return me.scales.x(d.mz);
    }).attr("y1", this.scales.y(0));
  },
  initScales: function() {
    this.ranges.x.min = 0;
    this.ranges.x.max = d3.max(this.data, function(r) { return r.mz; });
    this.ranges.y.min = 0;
    this.ranges.y.max = d3.max(this.data, function(r) { return r.intensity; });
    this.scales.x = d3.scale.linear().domain([this.ranges.x.min, this.ranges.x.max]).range([0, this.chartWidth]);
    this.scales.y = d3.scale.linear().domain([this.ranges.y.min, this.ranges.y.max]).range([this.chartHeight, 0]);
  },
  initAxes: function() {
    this.axes.x = d3.svg.axis().scale(this.scales.x).ticks(this.ticks.x);
    this.axes.y = d3.svg.axis().scale(this.scales.y).ticks(this.ticks.y).orient("left").tickFormat(d3.format('.2e'));
  },
  onDataReady: function() {
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
    this.svg.selectAll("line.mspeak")
    .data(this.data)
    .enter().append("svg:line")
    .attr("class", "mspeak")
    .attr("x1", function(d) { return me.scales.x(d.mz); })
    .attr("y2", function(d) { return me.scales.y(d.intensity); })
    .attr("y1", this.chartHeight)
    .attr("x2", function(d) { return me.scales.x(d.mz); })
    ;

    if (this.hasMarkers()) {
      this.onMarkersReady();
    }
  },
  setData: function(data) {
	  this.svg.selectAll('.axis').remove();
	  this.svg.selectAll('line.mspeak').remove();
    this.clearPeakSelection();
    this.svg.selectAll('.marker').remove();
    this.svg.selectAll('.emptytext').remove();
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
      me.fireEvent('selectpeak', mz);
      me.selectedpeak = mz;
    } else {
      me.fireEvent('unselectpeak', mz);
      me.selectedpeak = -1;
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
    this.markerSelect(function(d) {
      return (d.mz == mz);
    });
    this.selectedpeak = mz;
  },
  /**
   * Clears any selected peaks
   */
  clearPeakSelection: function() {
    this.markerSelect(false);
    this.selectedpeak = -1;
  },
  /**
   *
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
      .attr('class', 'marker lowermarker')
      .attr("transform", function(d) { return "translate(" + me.scales.x(d.mz) + "," + me.scales.y(0) + ")"; })
      .attr("d", d3.svg.symbol().type('triangle-up').size(36) )
      .on('click', markerClick)
      .append("svg:title")
        .text(markerTitle)
    ;

    // upper markers
    this.svg.selectAll("path.uppermarker")
    .data(function() {return me.markers;})
    .enter().append("svg:path")
      .attr('class', 'marker uppermarker')
      .attr("transform", function(d) { return "translate(" + me.scales.x(d.mz) + "," + me.scales.y(me.ranges.y.max) + ")"; })
      .attr("d", d3.svg.symbol().type('triangle-down').size(36) )
      .on('click', markerClick)
      .append("svg:title")
        .text(markerTitle)
    ;
  }
});