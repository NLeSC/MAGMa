/**
 * Ext wrapper around a D3 visualization.
 * @class Esc.d3.Abstract
 * @author Stefan Verhoeven
 *
 * @private
 */
Ext.define('Esc.d3.Abstract', {
  extend: 'Ext.Panel',
  initComponent: function() {
    var defConfig = {
        plain: true,
        border: false,
        /**
         * @property {Object} scales Scales
         * @property {d3.scale} scales.x X scale
         * @property {d3.scale} scales.y Y scale
         */
        scales: {
          x: null,
          y: null
        },
        /**
         * @property {Object} axes Axes
         * @property {d3.svg.axis} axes.x X axis
         * @property {d3.svg.axis} axes.y Y axis
         */
        axes: {
          x: null,
          y: null
        },
        /**
         * @cfg {Array} data array of objects.
         */
        data: [],
        emptyText: '',
        /**
         * @property {Object} ranges Range of axes
         * @property {Object} ranges.x Range of X axis
         * @property {Number} ranges.x.min Minimal x value. Default is 0.
         * @property {Number} ranges.x.min Maximal x value.
         * @property {Object} ranges.y Range of Y axis
         * @property {Number} ranges.y.min Minimal y value. Default is 0.
         * @property {Number} ranges.y.min Maximal y value.
         */
        ranges: {
          x: {
            min: 0,
            max: 0
          },
          y: {
            min: 0,
            max: 0
          },
        },
        /**
         * @cfg {Array} axesPadding Padding around axes. [top, right, left, bottom]
         */
        axesPadding: [10, 10, 38, 80],
        /**
         * @cfg {Object} ticks Number of ticks on axes.
         * @cfg {Number} [ticks.x=10] Number of ticks on x-axis.
         * @cfg {Number} [ticks.y=4] Number of ticks on y-axis.
         */
        ticks: {
          x:10,
          y:4
        },
        chartWidth: 0,
        chartHeight: 0
    };

    Ext.applyIf(this, defConfig);

    this.on('resize', this.onResize, this);
    // if width/height where unknown|tiny during onRender
    // they are known after layout so call initSvg again to
    // make sure svg is initialized
    this.on('afterlayout', function() {
      this.initSvg();
    }, this);
    this.callParent(arguments);
  },
  initSvg: function() {
    if (this.svg) {
      // dont reinit
      return;
    }
    var padding = this.axesPadding; // top right bottom left
    this.svg = d3.select(this.body.dom)
      .append('svg:svg')
      .attr('width',this.body.getWidth()).attr('height',this.body.getHeight())
      .attr('viewBox','0 0 '+this.body.getWidth()+' '+this.body.getHeight())
  //      .attr("preserveAspectRatio", "xMaxYMax meet")
      .attr("preserveAspectRatio", "none")
      .append('svg:g').attr('transform','translate('+padding[3]+','+padding[0]+')')
    ;
    // add
    this.svg.append('svg:rect')
        .attr('width',this.body.getWidth()).attr('height',this.body.getHeight())
        .attr('fill','none').attr('stroke', 'none')
    	.call(d3.behavior.zoom().on("zoom", this.redraw.bind(this) ))
    	.attr("pointer-events", "all")
    ;
    this.chartWidth = this.body.getWidth() - this.axesPadding[3] - this.axesPadding[1];
    this.chartHeight = this.body.getHeight() - this.axesPadding[0] - this.axesPadding[2];
    if (this.hasData()) {
      this.onDataReady();
    } else {
      this.onDataEmpty();
    }
  },
  hasData: function() {
    return (this.data.length>0);
  },
  onDataEmpty: function() {
    if (this.emptyText) {
       this.svg.append('svg:text')
         .attr('class', 'emptytext')
         .attr('x', this.chartWidth/2)
         .attr('y', this.chartHeight/2)
         .attr("dy", "1em")
         .attr("text-anchor", "middle")
         .attr('fill','gray')
         .text(this.emptyText);
    }
  },
  onResize: function(me, width, height) {
    // find svg tag and adjust w and h
    var s = d3.select(me.body.dom).select('svg');
    s.attr('width', width);
    s.attr('height', height);
  },
  onRender: function() {
    this.callParent(arguments);
    // width/height are very low when in this inside a % layout,
    // use afterlayout to get width/height
    if (this.body.getWidth() > 2 && this.body.getWidth() > 2) {
      this.initSvg();
    }
  },
  /**
   * Called after chart has been zoomed or panned.
   * See d3.behavior.zoom on zoom event.
   *
   * @method
   */
  redraw: Ext.emptyFn,
  /**
   * Prepare this.ranges and this.scales based on this.data .
   * Called by onDataReady.
   *
   * @method
   */
  initScales: Ext.emptyFn,
  /**
   * Prepare this.axes based on this.scales .
   * Called by onDataReady.
   *
   * @method
   */
  initAxes: Ext.emptyFn,
  /**
   * Plots this.data in this.svg.
   */
  onDataReady: function() {
    this.initScales();
    this.initAxes();
  },
  /**
   * Sets data and rerenders canvas
   * @param data
   */
  setData: function(data) {
    this.data = data;
    if (this.hasData()) {
      this.onDataReady();
    } else {
      this.onDataEmpty();
    }
  }
});