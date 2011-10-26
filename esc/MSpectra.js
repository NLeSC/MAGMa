/**
 * A Mass spectra viewer
 * @class Ext.esc.MSpectra
 * @extends Ext.Panel
 * @author Stefan Verhoeven
 */
Ext.define('Ext.esc.MSpectra', {
  extend: 'Ext.esc.AbstractD3',
  alias: 'widget.mspectra',
  /**
   * @cfg {Array} data array of objects with mz, intensity.
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
          .text('M/z')
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

    // basepeakintensity of each scan as vertical line
    this.svg.selectAll("line.mspeak")
    .data(this.data)
    .enter().append("svg:line")
    .attr("class", "mspeak")
    .attr("x1", function(d) { return me.scales.x(d.mz); })
    .attr("y2", function(d) { return me.scales.y(d.intensity); })
    .attr("y1", this.chartHeight)
    .attr("x2", function(d) { return me.scales.x(d.mz); })
    .style("stroke-width", '1px')
    .style("stroke", 'black')
    ;
  }
});