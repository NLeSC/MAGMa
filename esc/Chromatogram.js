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
        selectedScan: -1,
        /**
         * @cfg {String} selectedScanCls The CSS class applied to markers of a selected scan.
         */
        selectedScanCls: 'selectedscan',
        markers: [],
        chartWidth: 0,
        chartHeight: 0
    };

    Ext.applyIf(this, defConfig);

    this.callParent(arguments);

    this.addEvents(
      /**
       * @event selectscan
       * Fires when user clicks on scan marker (triangle)
       * @param {Int} scanid Scan identifier
       */
      'selectscan',
      /**
       * @event unselectscan
       * Fires when user clicks on a selected scan marker (triangle) which unselects it.
       * @param {Int} scanid Scan identifier
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
  },
  onDataReady: function() {
    var me = this;
    this.initScales();
    this.axes.x = d3.svg.axis().scale(this.scales.x).ticks(this.ticks.x);
    this.axes.y = d3.svg.axis().scale(this.scales.y).ticks(this.ticks.y).orient("left");

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
      .attr('stroke-dasharray','5,5')
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
    .style("stroke", function(d) {
      if (d.hashit) {
        return "#bbb";
      } else {
        return "#eee";
      }
    })
    ;

    // line drapped over peaks
    this.line = d3.svg.line()
    .interpolate('linear')
    .x(function(d) { return me.scales.x(d.rt); })
    .y(function(d) { return me.scales.y(d.intensity); });
    this.svg.append("svg:path")
    .attr("class", "line")
    .attr("d", this.line(this.data))
    .on('click', function(d) {
      console.log(d,d3.event);
    })
    ;
  },
  setData: function(data) {
	  this.svg.selectAll('.axis').remove();
	  this.svg.selectAll('.peak').remove();
	  this.svg.selectAll('.line').remove();
	  this.svg.selectAll('.'+this.cutoffCls).remove();
	  this.selectedscan = -1;
	  this.svg.selectAll('.marker').remove();
	  this.svg.selectAll('.emptytext').remove();
	  this.data = data;
    if (this.hasData()) {
      this.onDataReady();
    } else {
      this.onDataEmpty();
    }
  },
  onToggleMarker: function(scanid) {
    var me = this;
    this.markerSelect(function(e) {
      return (scanid == e.id && me.selectedscan != e.id);
    });
    if (scanid != me.selectedscan) {
      me.fireEvent('selectscan', scanid);
      me.selectedscan = scanid;
    } else {
      me.fireEvent('unselectscan', scanid);
      me.selectedscan = -1;
    }
  },
  /**
   * selects markers based on f returning true or false.
   * @param f function
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
  },
  /**
   * Clears any selected scans
   */
  clearScanSelection: function() {
    this.markerSelect(false);
    this.selectedscan = -1;
  },
  setMarkers: function(data) {
    this.clearScanSelection();
    this.svg.selectAll('.marker').remove();
    this.markers = data;
    this.onMarkersReady();
  },
  onMarkersReady: function() {
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
      .style("cursor", "pointer")
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
      .style("cursor", "pointer")
      .on('click', markerClick)
      .append("svg:title")
        .text(markerTitle)
    ;
  }
});