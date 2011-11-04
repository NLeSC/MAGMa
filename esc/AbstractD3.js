/**
 * Ext wrapper around a D3 visualization.
 * @class Ext.esc.AbstractD3
 * @extends Ext.Panel
 * @author Stefan Verhoeven
 */
Ext.define('Ext.esc.AbstractD3', {
  extend: 'Ext.Panel',

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
         * @cfg {Array} data array of objects.
         */
        data: [],
        emptyText: '',
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
        axesPadding: [16, 5, 38, 80],
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
        chartWidth: 0,
        chartHeight: 0
    };

    Ext.applyIf(this, defConfig);

    this.callParent(arguments);
  },
  initSvg: function() {
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
  onRender: function() {
    this.callParent(arguments);
    // width/height are very low when in this inside a % layout,
    // use afterlayout to get width/height
    if (this.body.getWidth() > 2 && this.body.getWidth() > 2) {
      this.initSvg();
    }
    this.on('afterlayout', function() {
      if (!this.svg) {
        this.initSvg();
      }
    });
    this.on('resize', function(t,width, height) {
      // find svg tag and adjust w and h
      var s = d3.select(t.body.dom).select('svg');
      s.attr('width', this.body.getWidth());
      s.attr('height', this.body.getHeight());
    });
  },
  redraw: Ext.emptyFn,
  initScales: Ext.emptyFn,
  onDataReady: Ext.emptyFn,
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