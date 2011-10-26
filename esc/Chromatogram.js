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
          x: null,
          y: null
        },
        axes: {
          x: null,
          y: null
        },
        data: [],
        ranges: {
          x: { min: 0, max: 0},
          y: { min: 0, max: 0},
        },
        axesPadding: [16, 5, 50, 80],
        cutoff: 2000000,
        ticks: { x:10, y:4 },
        selectedScan: -1,
        chartWidth: 0,
        chartHeight: 0
    };

    Ext.applyIf(this, defConfig);

    this.callParent(arguments);
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
      .attr('class','cutoffline')
      .attr('x1',0)
      .attr('x2',this.chartWidth)
      .attr('y1',this.scales.y(this.cutoff))
      .attr('y2',this.scales.y(this.cutoff))
      .attr('stroke-dasharray','5,5')
    ;

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

    var peaks_with_hits = this.data.filter(function(d) {
      return d.hashit;
    });

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
  markerSelect: function(f) {
    this.svg.selectAll("path.uppermarker")
    .classed("selectedscan", f)
    ;
    this.svg.selectAll("path.lowermarker")
      .classed("selectedscan", f)
    ;
  },
  selectScans: function(scanids) {
    this.markerSelect(function(d) {
      return (d.id in scanids);
    });
  },
  clearScanSelection: function(scanids) {
    this.markerSelect(false);
    this.selectedscan = -1;
  },
  resetZoom: function() {
    // reset d3.behavior.zoom, by overwriting
    d3.select(this.body.dom).select('svg').call(
        d3.behavior.zoom().on("zoom", this.redraw.bind(this) )
    );
    this.initScales();
    this.redraw();
  },
  setData: function(data) {
    this.data = data;
    this.onDataReady();
  }
});