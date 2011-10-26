/**
 * @class Ext.esc.Chromatogram
 * @extends Ext.Panel
 * @author Stefan Verhoeven
 */
Ext.define('Ext.esc.Chromatogram', {
  extend: 'Ext.Panel',
  alias: 'widget.chromatogram',

  initComponent: function() {
    var defConfig = {
        plain: true,
        border: false,
        scales: {
          /**
           * @cfg {d3.scale} scales.x X scale
           */
          x: null,
          /**
           * @cfg {d3.scale} scales.y Y scale
           */
          y: null
        },
        axes: {
          /**
           * @cfg {d3.svg.axis} axes.x X axis
           */
          x: null,
          /**
           * @cfg {d3.svg.axis} axes.y Y axis
           */
          y: null
        },
        /**
         * @cfg {Array} data array of objects with id, rt, intensity and hashit properties.
         */
        data: [],
        ranges: {
          x: {
            /**
             * @cfg {Number} ranges.x.min Minimal x value. Default is 0.
             */
            min: 0,
            /**
             * @cfg {Number} ranges.x.min Maximal x value. Default is max rt in this.data.
             */
            max: 0
          },
          y: {
            /**
             * @cfg {Number} ranges.x.min Minimal y value. Default is 0.
             */
            min: 0,
            /**
             * @cfg {Number} ranges.x.min Maximal y value. Default is max intensity in this.data.
             */
            max: 0
          },
        },
        /**
         * @cfg {Array} axesPadding Padding around axes. [top, right, left, bottom]
         */
        axesPadding: [16, 5, 50, 80],
        /**
         * @cfg {Number} cutoff Intensity under which scans where disregarded
         */
        cutoff: 2000000,
        /**
         * @cfg {String} cutoffCls The CSS class applied to cutoff line.
         */
        cutoffCls: 'cutoffline',
        ticks: {
          /**
           * @cfg {Number} ticks.x=10 Number of ticks on x-axis.
           */
          x:10,
          /**
           * @cfg {Number} ticks.y=10 Number of ticks on y-axis.
           */
          y:4
        },
        selectedScan: -1,
        /**
         * @cfg {String} selectedScanCls The CSS class applied to markers of a selected scan.
         */
        selectedScanCls: 'selectedscan',
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
  onRender: function() {
    this.callParent(arguments);
    var padding = this.axesPadding; // top right bottom left
    this.svg = d3.select(this.body.dom)
      .append('svg:svg')
      .attr('width',this.body.getWidth()).attr('height',this.body.getHeight())
      .attr('viewBox','0 0 '+this.body.getWidth()+' '+this.body.getHeight())
//      .attr("preserveAspectRatio", "xMaxYMax meet")
      .attr("preserveAspectRatio", "none")
      .attr("pointer-events", "all")
      .call(d3.behavior.zoom().on("zoom", this.redraw.bind(this) ))
      .append('svg:g').attr('transform','translate('+padding[3]+','+padding[0]+')')
      ;
    this.chartWidth = this.body.getWidth() - this.axesPadding[3] - this.axesPadding[1];
    this.chartHeight = this.body.getHeight() - this.axesPadding[0] - this.axesPadding[2];
    if (this.data.length) {
      this.onDataReady();
    }
    this.on('resize', function(t,width, height) {
      // find svg tag and adjust w and h
      var s = d3.select(t.body.dom).select('svg');
      s.attr('width', this.body.getWidth());
      s.attr('height', this.body.getHeight());
    });
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
    this.svg.selectAll("path.lowermarker")
      .attr("transform", function(d) { return "translate(" + me.scales.x(d.rt) + "," + me.scales.y(0) + ")"; });
    this.svg.selectAll("path.uppermarker")
     .attr("transform", function(d) { return "translate(" + me.scales.x(d.rt) + "," + me.scales.y(me.ranges.y.max) + ")"; });
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

    // add markers to peaks which have hit
    var peaks_with_hits = this.data.filter(function(d) {
      return d.hashit;
    });

    // lower markers
    this.svg.selectAll("path.lowermarker")
    .data(peaks_with_hits)
    .enter().append("svg:path")
      .attr('class', 'marker lowermarker')
      .attr("transform", function(d) { return "translate(" + me.scales.x(d.rt) + "," + me.scales.y(0) + ")"; })
      .attr("d", d3.svg.symbol().type('triangle-up').size(36) )
      .style("cursor", "pointer")
      .on('click', function(d) {
        me.onToggleMarker(d.id);
      })
      .append("svg:title")
        .text(function(d) { return 'Scan#'+d.id; })
    ;

    // upper markers
    this.svg.selectAll("path.uppermarker")
    .data(peaks_with_hits)
    .enter().append("svg:path")
      .attr('class', 'marker uppermarker')
      .attr("transform", function(d) { return "translate(" + me.scales.x(d.rt) + "," + me.scales.y(me.ranges.y.max) + ")"; })
      .attr("d", d3.svg.symbol().type('triangle-down').size(36) )
      .style("cursor", "pointer")
      .on('click', function(d) {
        me.onToggleMarker(d.id);
      })
      .append("svg:title")
        .text(function(d) { return 'Scan#'+d.id; })
    ;
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
      return (d.id in scanids);
    });
  },
  /**
   * Clears any selected scans
   */
  clearScanSelection: function() {
    this.markerSelect(false);
    this.selectedscan = -1;
  },
  /**
   * Reset zoom.
   */
  resetZoom: function() {
    // reset d3.behavior.zoom, by overwriting
    d3.select(this.body.dom).select('svg').call(
        d3.behavior.zoom().on("zoom", this.redraw.bind(this) )
    );
    this.initScales();
    this.redraw();
  },
  /**
   * Stuffs data into chromatogram.
   * @param data
   */
  setData: function(data) {
    this.data = data;
    this.onDataReady();
  }
});