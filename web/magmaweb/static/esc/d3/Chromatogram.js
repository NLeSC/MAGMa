/**
 * A chromatogram viewer
 * @class Esc.d3.Chromatogram
 * @extends Esc.d3.Abstract
 * @author Stefan Verhoeven
 *
 * A example with 2 peaks, a marker at the first peak, an extracted ion chromatogram and a select of a scan:
 *
 *     @example
 *     var chromatogram = Ext.create('Esc.d3.Chromatogram', {
 *       renderTo: Ext.getBody(),
 *       title: 'Chromatogram',
 *       width: 400, height: 300,
 *       axesPadding: [16, 5, 58, 80],
 *       data: [{rt:1, intensity: 100, id:4}, {rt:2, intensity: 50, id:5}],
 *       markers: [{rt:1, intensity: 100, id:4}],
 *       cutoff: 10
 *     });
 *     chromatogram.setExtractedIonChromatogram([{rt:1, intensity: 45, id:4}, {rt:2, intensity: 60, id:5}]);
 *     chromatogram.selectScan(4);
 *
 * Note! This example requires d3.v3.js to be sourced.
 */
Ext.define('Esc.d3.Chromatogram', {
  extend: 'Esc.d3.Abstract',
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
       * When no scans are selected then it is set to false.
       * @property {Number}
       * @readonly
       */
      selectedScan: false,
      /**
       * @cfg {String} selectedScanCls The CSS class applied to markers of a selected scan.
       */
      selectedScanCls: 'selected',
      // array of {rt:,intensity:, id:}
      markers: [],
      chartWidth: 0,
      chartHeight: 0,
      // array of {rt:,intensity:}
      moleculedata: []
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
      'unselectscan',
      /**
       * @event
       * Fires when mouse is moved over a vertical line of the scan
       * @param {Object} scan
       * @param {Number} scan.id Scan identifier
       * @param {Number} scan.rt Retention time
       * @param {Number} scan.intensity Intensity of base peak.
       * @param {Number} scan.moleculeintensity If extracted ion chromatogram is set then returns intensity of molecule in scan.
       */
      'mouseoverscan'
    );
  },
  /**
   * @inheritdoc Esc.d3.Abstract#onZoom
   */
  onZoom: function() {
    this.callParent(arguments);
    var me = this;

    this.svg.select("path.line").attr('d', this.line(this.data));
    this.svg.select("path.moleculeline").attr('d', this.line(this.moleculedata));
    if (this.markers.length) {
      this.svg.selectAll("path.lowermarker")
        .attr("transform", function(d) {
          return "translate(" + me.scales.x(d.rt) + "," + (me.scales.y(0) + 4) + ")";
        });
      this.svg.selectAll("path.uppermarker")
        .attr("transform", function(d) {
          return "translate(" + me.scales.x(d.rt) + "," + (me.scales.y(me.ranges.y.max) - 4) + ")";
        });
    }
    this.svg.selectAll("line.peak")
      .attr("x1", function(d) {
        return me.scales.x(d.rt);
      })
      .attr("y2", function(d) {
        return me.scales.y(d.intensity);
      })
      .attr("y1", function(d) {
        return me.scales.y(0);
      })
      .attr("x2", function(d) {
        return me.scales.x(d.rt);
      });

    this.svg.selectAll("line." + this.cutoffCls)
      .attr('x1', 0)
      .attr('x2', this.chartWidth)
      .attr('y1', this.scales.y(this.cutoff))
      .attr('y2', this.scales.y(this.cutoff));
  },
  initScales: function() {
    this.callParent(arguments);
    this.ranges.x.min = 0;
    this.ranges.x.max = d3.max(this.data, function(r) {
      return r.rt;
    });
    this.ranges.y.min = 0;
    this.ranges.y.max = d3.max(this.data, function(r) {
      return r.intensity;
    });
    this.scales.x = d3.scale.linear().domain([this.ranges.x.min, this.ranges.x.max]).range([0, this.chartWidth]);
    this.scales.y = d3.scale.linear().domain([this.ranges.y.min, this.ranges.y.max]).range([this.chartHeight, 0]);
  },
  initAxes: function() {
    this.callParent(arguments);
    var me = this;
    /**
     * @property {d3.svg.line} line Line factory for basepeakintensity and extractedionchromatogram
     */
    this.line = d3.svg.line()
      .interpolate('linear')
      .x(function(d) {
        return me.scales.x(d.rt);
      })
      .y(function(d) {
        return me.scales.y(d.intensity);
      });

    var nrxticks = this.ticks.x;
    if (this.chartWidth < 25 * (6 + 2)) {
      nrxticks = 2;
    }
    if (this.chartWidth < 25 * 2) {
      nrxticks = 0;
    }
    this.axes.x = d3.svg.axis().scale(this.scales.x).ticks(nrxticks);

    var nryticks = this.ticks.y;
    if (this.chartHeight < 16 * (1 + 2)) {
      nryticks = 2;
    }
    if (this.chartHeight < 16 * 2) {
      nryticks = 0;
    }
    this.axes.y = d3.svg.axis().scale(this.scales.y).ticks(nryticks).orient("left").tickFormat(d3.format('.2e'));
  },
  draw: function() {
    this.callParent(arguments);
    var me = this;

    // Add the x-axis.
    this.svg.append("svg:g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + this.chartHeight + ")")
      .call(this.axes.x)
      .append("svg:text")
      .attr("x", this.chartWidth / 2).attr("y", 30)
      .attr("text-anchor", "middle")
      .text('Retention time (min)');

    // Add the y-axis.
    this.svg.append("svg:g")
      .attr("class", "y axis")
      .call(this.axes.y)
      .append("svg:text")
      .attr("y", this.chartHeight / 2)
      .attr("x", -5)
      .attr("text-anchor", "middle")
      .attr("transform", "rotate(-90," + -5 + "," + this.chartHeight / 2 + ")")
      .text('Intensity');

    // cutoff
    this.svg.append("svg:line")
      .attr('class', this.cutoffCls)
      .attr('x1', 0)
      .attr('x2', this.chartWidth)
      .attr('y1', this.scales.y(this.cutoff))
      .attr('y2', this.scales.y(this.cutoff));

    // basepeakintensity of each scan as vertical line
    this.svg.selectAll("line.peak")
      .data(this.data)
      .enter().append("svg:line")
      .attr("class", "peak")
      .classed('assigned', function(d) {
        return d.ap > 0;
      })
      .attr("x1", function(d) {
        return me.scales.x(d.rt);
      })
      .attr("y2", function(d) {
        return me.scales.y(d.intensity);
      })
      .attr("y1", function(d) {
        return me.scales.y(0);
      })
      .attr("x2", function(d) {
        return me.scales.x(d.rt);
      })
      .on('mouseover', function(scan) {
        // fetch intensity of molecule if available
        if (me.moleculedata.length) {
          scan.moleculeintensity = me.moleculedata.filter(function(d) {
            return (scan.rt == d.rt);
          })[0].intensity;
        }
        me.fireEvent('mouseoverscan', scan);
      });

    // line drapped over peaks
    this.svg.append("svg:path")
      .attr("class", "line")
      .attr("d", this.line(this.data));

    if (this.hasMarkers()) {
      this.onMarkersReady();
    }
  },
  undraw: function() {
    this.svg.selectAll('.axis').remove();
    this.svg.selectAll('.peak').remove();
    this.svg.selectAll('.line').remove();
    this.svg.selectAll('.' + this.cutoffCls).remove();
    this.clearScanSelection();
    this.svg.selectAll('.marker').remove();
    this.moleculedata = [];
    this.svg.selectAll('path.moleculeline').remove();
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
      if (me.selectedScan) {
        me.fireEvent('unselectscan', me.selectedScan);
      }
      me.fireEvent('selectscan', scanid);
      me.selectedScan = scanid;
    } else {
      me.fireEvent('unselectscan', scanid);
      me.selectedScan = false;
    }
  },
  /**
   * selects markers based on f returning true or false.
   * @param f function
   * @private
   */
  markerSelect: function(f) {
    this.svg.selectAll("path.uppermarker")
      .classed(this.selectedScanCls, f);
    this.svg.selectAll("path.lowermarker")
      .classed(this.selectedScanCls, f);
  },
  /**
   * Select scan by their id
   * When there is another scan already selected the unselectscan event is fired.
   * @param {Number} scanid Scan identifier.
   * @param {Boolean} [silent=false] Passing true will supress the 'unselectscan' event from being fired.
   */
  selectScan: function(scanid, silent) {
    if (this.selectedScan && this.selectedScan != scanid && silent !== true) {
      this.fireEvent('unselectscan', this.selectedScan);
    }
    this.markerSelect(function(d) {
      return (scanid == d.id);
    });
    if (scanid) {
      this.selectedScan = scanid;
    } else {
      this.selectedScan = false;
    }
  },
  /**
   * Clears any selected scans
   */
  clearScanSelection: function() {
    this.markerSelect(false);
    this.selectedScan = false;
  },
  /**
   * Set markers on rt's which can be selected.
   * @param {Array} data array of markers.
   *
   * When a scan has been selected scan it is reselected if scan is still marked
   */
  setMarkers: function(data) {
    var selectedScan;
    if (this.selectedScan) {
      selectedScan = this.selectedScan;
    }
    this.clearScanSelection();
    this.svg.selectAll('.marker').remove();
    this.markers = data;
    this.onMarkersReady();
    if (selectedScan) {
      this.selectScan(selectedScan);
    }
  },
  hasMarkers: function() {
    return (this.markers.length > 0);
  },
  onMarkersReady: function() {
    // can not add markers if there is no data
    if (!this.hasData()) {
      return;
    }

    var me = this;

    function markerTitle(d) {
      return 'Scan#' + d.id;
    }

    function markerClick(d) {
      me.onToggleMarker(d.id);
    }

    // lower markers
    this.svg.selectAll("path.lowermarker")
      .data(function() {
        return me.markers;
      })
      .enter().append("svg:path")
      .attr('class', 'marker lowermarker annotated')
      .attr("transform", function(d) {
        return "translate(" + me.scales.x(d.rt) + "," + (me.scales.y(0) + 4) + ")";
      })
      .attr("d", d3.svg.symbol().type('triangle-up').size(32))
      .on('click', markerClick)
      .append("svg:title")
      .text(markerTitle);

    // upper markers
    this.svg.selectAll("path.uppermarker")
      .data(function() {
        return me.markers;
      })
      .enter().append("svg:path")
      .attr('class', 'marker uppermarker annotated')
      .attr("transform", function(d) {
        return "translate(" + me.scales.x(d.rt) + "," + (me.scales.y(me.ranges.y.max) - 4) + ")";
      })
      .attr("d", d3.svg.symbol().type('triangle-down').size(32))
      .on('click', markerClick)
      .append("svg:title")
      .text(markerTitle);
  },
  /**
   * Overlay the extracted ion chromatogram of a molecule on the chromatogram.
   * @param data Array of rt and max intensity of a molecule
   */
  setExtractedIonChromatogram: function(data) {
    this.moleculedata = data;
    this.svg.selectAll('path.moleculeline').remove();
    this.svg.append('svg:path')
      .attr('class', 'moleculeline')
      .attr('d', this.line(this.moleculedata));
  }
});
