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
        axesPadding: [10, 10, 38, 80],
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
   * Adds svg canvas to component.
   */
  initSvg: function() {
    if (this.svg) {
      // dont reinit
      return;
    }

    var m = this.axesPadding;
    this.svg = d3.select('#'+this.id).append('svg:svg').append("g")
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
    this.svg.select(".x.axis").call(this.axes.x);
  },
  /**
   * Prepare {@link #ranges}, {@link #scales} and {@link #axes} based on {@link #data} .
   * @template
   */
  initScales: function() {
    this.chartWidth = this.getWidth() - this.axesPadding[1] - this.axesPadding[3];
    this.chartHeight = this.getHeight() - this.axesPadding[0] - this.axesPadding[2];
  },
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
  /**
   * @protected
   * Applies zoom behavior on x-axis to canvas
   */
  initZoom: function() {
    // update zoomer
    this.svg.select('rect.zoomer').call(
        d3.behavior.zoom().x(this.scales.x).on("zoom", this.onZoom.bind(this))
    );
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