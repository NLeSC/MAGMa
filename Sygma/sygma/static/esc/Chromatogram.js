/**
 * A chromatogram viewer
 * @class Ext.esc.Chromatogram
 * @extends Ext.Panel
 * @author Stefan Verhoeven
 */
Ext.define('Ext.esc.Chromatogram', {
  extend: 'Ext.esc.AbstractD3',
  alias: 'widget.chromatogram',

  initComponent: function() {
    var defConfig = {
        /**
         * @cfg {Array} data array of objects with scanid, rt and intensity properties.
         */

        /**
         * @cfg {Number} cutoff Intensity under which scans where disregarded
         */
        cutoff: 2000000,
        /**
         * @cfg {String} cutoffCls The CSS class applied to cutoff line.
         */
        cutoffCls: 'cutoffline',
        /**
         * Scan identifier of selected scan.
         * When no scans are selected then it is set to -1.
         * @prop {Number}
         * @readonly
         */
        selectedScan: -1,
        /**
         * @cfg {String} selectedScanCls The CSS class applied to markers of a selected scan.
         */
        selectedScanCls: 'selectedscan',
        // array of {rt:,intensity:, id:}
        markers: [],
        chartWidth: 0,
        chartHeight: 0,
        // array of {rt:,intensity:}
        metabolitedata: [],
    };

    Ext.applyIf(this, defConfig);

    this.callParent(arguments);

    this.addEvents(
      /**
       * @event selectscan
       * Fires when user clicks on scan marker (triangle)
       * @param {Number} scanid Scan identifier
       */
      'selectscan',
      /**
       * @event unselectscan
       * Fires when user clicks on a selected scan marker (triangle) which unselects it.
       * @param {Number} scanid Scan identifier
       */
      'unselectscan'
    );
  },
  redraw: function() {
    var me = this;
    if (d3.event && d3.event.translate[0] != 0 && d3.event.translate[1] != 0) {
      // pan and zoom x axis
      d3.event.transform(this.scales.x);
    }
    this.svg.select(".x.axis").call(this.axes.x);
    // do not scale y axis
    //svg.select(".y.axis").call(yAxis);
    this.svg.select("path.line").attr('d', this.line(this.data));
    this.svg.select("path.metaboliteline").attr('d', this.line(this.metabolitedata));
    if (this.markers.length) {
      this.svg.selectAll("path.lowermarker")
        .attr("transform", function(d) { return "translate(" + me.scales.x(d.rt) + "," + me.scales.y(0) + ")"; });
      this.svg.selectAll("path.uppermarker")
       .attr("transform", function(d) { return "translate(" + me.scales.x(d.rt) + "," + me.scales.y(me.ranges.y.max) + ")"; });
    }
    this.svg.selectAll("line.peak").attr("y2", function(d) {
      return me.scales.y(d.intensity);
    }).attr("x1", function(d) {
      return me.scales.x(d.rt);
    }).attr("x2", function(d) {
      return me.scales.x(d.rt);
    }).attr("y1", this.scales.y(0));
  },
  initScales: function() {
    this.ranges.x.min = 0;
    this.ranges.x.max = d3.max(this.data, function(r) { return r.rt; });
    this.ranges.y.min = 0;
    this.ranges.y.max = d3.max(this.data, function(r) { return r.intensity; });
    this.scales.x = d3.scale.linear().domain([this.ranges.x.min, this.ranges.x.max]).range([0, this.chartWidth]);
    this.scales.y = d3.scale.linear().domain([this.ranges.y.min, this.ranges.y.max]).range([this.chartHeight, 0]);
    var me = this;
    /**
     * Line factory for basepeakintensity and extractedionchromatogram
     * @cfg line
     */
    this.line = d3.svg.line()
    .interpolate('linear')
    .x(function(d) { return me.scales.x(d.rt); })
    .y(function(d) { return me.scales.y(d.intensity); });
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
          .text('Retention time (s)')
    ;

    // Add the y-axis.
    this.svg.append("svg:g")
        .attr("class", "y axis")
        .call(this.axes.y)
        .append("svg:text")
          .attr("y",this.chartHeight/2)
          .attr("x",-5)
          .attr("text-anchor", "middle")
          .attr("transform", "rotate(-90,"+-5+","+this.chartHeight/2+")" )
          .text('Intensity')
    ;

    // cutoff
    this.svg.append("svg:line")
      .attr('class', this.cutoffCls)
      .attr('x1',0)
      .attr('x2',this.chartWidth)
      .attr('y1',this.scales.y(this.cutoff))
      .attr('y2',this.scales.y(this.cutoff))
    ;

    // basepeakintensity of each scan as vertical line
    this.svg.selectAll("line.peak")
    .data(this.data)
    .enter().append("svg:line")
    .attr("class", "peak")
    .attr("x1", function(d) { return me.scales.x(d.rt); })
    .attr("y2", function(d) { return me.scales.y(d.intensity); })
    .attr("y1", this.chartHeight)
    .attr("x2", function(d) { return me.scales.x(d.rt); })
    ;

    // line drapped over peaks
    this.svg.append("svg:path")
    .attr("class", "line")
    .attr("d", this.line(this.data))
    ;

    if (this.hasMarkers()) {
      this.onMarkersReady();
    }
  },
  setData: function(data) {
	  this.svg.selectAll('.axis').remove();
	  this.svg.selectAll('.peak').remove();
	  this.svg.selectAll('.line').remove();
	  this.svg.selectAll('.'+this.cutoffCls).remove();
	  this.clearScanSelection();
	  this.svg.selectAll('.marker').remove();
	  this.svg.selectAll('.emptytext').remove();
	  this.metabolitedata = [];
    this.svg.selectAll('path.metaboliteline').remove();
    this.callParent(arguments);
  },
  /**
   * @private
   */
  onToggleMarker: function(scanid) {
    var me = this;
    this.markerSelect(function(e) {
      return (scanid == e.id && me.selectedScan != e.id);
    });
    if (scanid != me.selectedScan) {
      me.fireEvent('selectscan', scanid);
      me.selectedScan = scanid;
    } else {
      me.fireEvent('unselectscan', scanid);
      me.selectedScan = -1;
    }
  },
  /**
   * selects markers based on f returning true or false.
   * @param f function
   * @private
   */
  markerSelect: function(f) {
    this.svg.selectAll("path.uppermarker")
      .classed(this.selectedScanCls, f)
    ;
    this.svg.selectAll("path.lowermarker")
      .classed(this.selectedScanCls, f)
    ;
  },
  /**
   * Select scans by their id
   * @param scanids Array of scan ids
   */
  selectScans: function(scanids) {
    this.markerSelect(function(d) {
      return !(scanids.indexOf(d.id) == -1);
    });
    if (scanids.length == 1) {
       this.selectedScan = scanids[0];
     } else if (scanids.length == 0) {
       this.selectedScan = -1;
     }
  },
  /**
   * Clears any selected scans
   */
  clearScanSelection: function() {
    this.markerSelect(false);
    this.selectedScan = -1;
  },
  /**
   * @param data array of markers.
   *
   * When a scan has been selected scan it is reselected if scan is still marked
   */
  setMarkers: function(data) {
  	if (this.selectedScan != -1) {
  		var selectedScan = this.selectedScan;
  	}
    this.clearScanSelection();
    this.svg.selectAll('.marker').remove();
    this.markers = data;
    this.onMarkersReady();
    if (selectedScan) {
    	this.selectScans([selectedScan]);
    }
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
      return 'Scan#'+d.id;
    }
    function markerClick(d) {
      me.onToggleMarker(d.id);
    }

    // lower markers
    this.svg.selectAll("path.lowermarker")
    .data(function() {return me.markers;})
    .enter().append("svg:path")
      .attr('class', 'marker lowermarker')
      .attr("transform", function(d) { return "translate(" + me.scales.x(d.rt) + "," + me.scales.y(0) + ")"; })
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
      .attr("transform", function(d) { return "translate(" + me.scales.x(d.rt) + "," + me.scales.y(me.ranges.y.max) + ")"; })
      .attr("d", d3.svg.symbol().type('triangle-down').size(36) )
      .on('click', markerClick)
      .append("svg:title")
        .text(markerTitle)
    ;
  },
  /**
   * @param data Array of rt and max intensity of a metabolite
   */
  setExtractedIonChromatogram: function(data) {
    this.metabolitedata = data;
    this.svg.selectAll('path.metaboliteline').remove();
    this.svg.append('svg:path')
      .attr('class','metaboliteline')
      .attr('d', this.line(this.metabolitedata) )
    ;
  }
});