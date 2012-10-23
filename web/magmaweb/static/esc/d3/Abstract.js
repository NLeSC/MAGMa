/**
 * Ext wrapper around a D3 visualization.
 * @class Esc.d3.Abstract
 * @author Stefan Verhoeven
 *
 * This is an abstract superclass and should not be used directly.
 * @private
 */
Ext.define('Esc.d3.Abstract', {
  extend: 'Ext.Component',
  initComponent: function() {
    var defConfig = {
        plain: true,
        /**
         * @property {Object} scales Scales used by axes
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
         * @property {Object} zoom Zoom on axes
         * @property {Boolean} zoom.x Enable zooming/panning on X axis. Default true.
         * @property {Boolean} zoom.y Enable zooming/panning on Y axis. Default false.
         */
        zoom: {
          x: true,
          y: false
        },
        /**
         * @cfg {Array} data array of objects.
         */
        data: [],
        /**
         * @cfg {String} emptyText The text to display in the view when there is no data to display.
         */
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
          }
        },
        /**
         * @cfg {Array} axesPadding Padding around axes. [top, right, left, bottom]
         */
        axesPadding: [10, 14, 38, 80],
        /**
         * @cfg {Object} ticks Number of ticks on axes.
         * @cfg {Number} [ticks.x=10] Number of ticks on x-axis.
         * @cfg {Number} [ticks.y=4] Number of ticks on y-axis.
         */
        ticks: {
          x: 10,
          y: 4
        },
        /**
         * @property {Number} chartWidth Width of area in canvas where data is plotted.
         * Is width of component minus {@link #axesPadding}.
         */
        chartWidth: 0,
        /**
         * @property {Number} chartHeight Height of area in canvas where data is plotted.
         * Is height of component minus {@link #axesPadding}.
         */
        chartHeight: 0
    };

    Ext.applyIf(this, defConfig);

    this.on('resize', this.onResize, this);
    // if width/height where unknown|tiny during onRender
    // they are known after layout so call initSvg again to
    // make sure svg is initialized
    this.callParent(arguments);
  },
  /**
   * @protected
   *
   * Depending if data is set will draw data on canvas or display emptytext.
   */
  onResize: function() {
    if (this.hasData()) {
      this.onDataReady();
    } else {
      this.onDataEmpty();
    }
  },
  /**
   * @protected
   *
   * Adds svg canvas to component.
   */
  initSvg: function() {
    if (this.svg) {
      // dont reinit
      return;
    }

    var m = this.axesPadding;
    this.svg = d3.select('#'+this.id).append('svg:svg')
        .attr('width', '100%').attr('height', '100%')
        .append("g")
        .attr("transform", "translate(" + m[3] + "," + m[0] + ")")
    ;

    // zoomer rect which captures mouse drags and mouse wheel events
    this.svg.append('svg:rect')
        .attr('class', 'zoomer')
        .attr('fill','none').attr('stroke', 'none')
        .attr('width', '100%')
        .attr('height', '100%')
        .attr("pointer-events", "all")
    ;

    this.onResize();
  },
  /**
   * Whether data is associated.
   * @return {Boolean}
   */
  hasData: function() {
    return (this.data.length>0);
  },
  /**
   * @protected
   * Clears canvas and prints emptytext if available
   */
  onDataEmpty: function() {
    this.undraw();
    if (this.emptyText) {
       this.svg.append('svg:text')
         .attr('class', 'emptytext')
         .attr('x', 0)
         .attr('y', 0)
         .attr("dx", "1em")
         .attr("dy", "1em")
         .attr("text-anchor", "begin")
         .attr('fill','gray')
         .text(this.emptyText);
    }
  },
  /**
   * Adds svg canvas to component dom
   */
  onRender: function() {
    this.callParent(arguments);
    this.initSvg();
  },
  /**
   * Called after chart has been zoomed or panned.
   * {@link #scales} have been updated.
   * Transform chart to fit new scales.
   *
   * @template
   */
  onZoom: function() {
  	if (d3.event) {
    	var translate = d3.event.translate;
    	var scale = d3.event.scale;
  	} else {
  		var translate = [0,0];
  		var scale = 1;
  	}

    if (this.zoom.x) {
      var x1 = this.scales.x;
      var x0 = this.scales0.x;
      x1.domain(x0.range().map(function(x) { return (x - translate[0]) / scale; }).map(x0.invert));
      this.svg.select(".x.axis").call(this.axes.x);
    }

    if (this.zoom.y) {
      var y1 = this.scales.y;
      var y0 = this.scales0.y;
      // ignore translate, only scale
      var r = y0.range().map(function(y) { return (y / scale ); });
      // r[1] is fixed, we want r[0] fixed so swap them
      r = [y0.range()[0], y0.range()[0]-r[0]];
      y1.domain(r.map(y0.invert));
      this.svg.select(".y.axis").call(this.axes.y);
    }
  },
  /**
   * Prepare {@link #ranges}, {@link #scales} based on {@link #data} .
   * @template
   */
  initScales: function() {
    this.chartWidth = this.getWidth() - this.axesPadding[1] - this.axesPadding[3];
    this.chartHeight = this.getHeight() - this.axesPadding[0] - this.axesPadding[2];
  },
  /**
   * Prepare {@link #axes} based on {@link #data} .
   */
  initAxes: Ext.emptyFn,
  /**
   * Render axis and data to canvas.
   * @template
   */
  draw: function() {
    this.initZoom();
  },
  /**
   * Remove rendered axis and data from canvas.
   * @template
   *
   */
  undraw: function() {
    this.svg.selectAll('text.emptytext').remove();
  },
  /**
   * @protected
   * Plots this.data in this.svg.
   */
  onDataReady: function() {
    this.initScales();
    this.initAxes();
    this.undraw();
    this.draw();
  },
  /**
   * Sets data and rerenders canvas
   * @param data Array
   */
  setData: function(data) {
    this.data = data;
    this.onResize();
  },
  zoomBehavior: function() {
    var zoom = d3.behavior.zoom().on("zoom", this.onZoom.bind(this));

    // store initial state of scales so zoom/translate can use it as reference
    this.scales0 = {};
    this.scales0.x = this.scales.x.copy();
    this.scales0.y = this.scales.y.copy();

    return zoom;
  },
  /**
   * @protected
   * Applies zoom behavior on x-axis to canvas
   */
  initZoom: function() {
    // update zoomer
    this.svg.select('rect.zoomer').call(this.zoomBehavior());
  },
  /**
   * Enable zooming/panning on an axis.
   * @param {String} axis Name of an axis, can be 'x' or 'y'.
   * @param {Boolean} enabled True to enable.
   */
  setZoom: function(axis, enabled) {
    this.zoom[axis] = enabled;
    this.initZoom();
  },
  /**
   * Resets scales back to their original ranges
   */
  resetScales: function() {
    this.scales.x.domain([this.ranges.x.min, this.ranges.x.max]);
    this.scales.y.domain([this.ranges.y.min, this.ranges.y.max]);
    this.initZoom();
    this.onZoom();
  }
});
