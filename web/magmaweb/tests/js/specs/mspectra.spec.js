describe('Esc.d3.MSpectra', function() {
  var data = [{
    mz: 1,
    intensity: 2
  }, {
    mz: 3,
    intensity: 4
  }];

  function mockSvg() {
    var svg = {
      selectAll: function() {
        return this;
      },
      select: function() {
        return this;
      },
      append: function() {
        return this;
      },
      attr: function() {
        return this;
      },
      call: function() {
        return this;
      },
      data: function() {
        return this;
      },
      enter: function() {
        return this;
      },
      remove: function() {
        return this;
      },
      classed: function() {
        return this;
      },
      on: function() {
        return this;
      },
      text: function() {
        return this;
      }
    };
    spyOn(svg, 'selectAll').andCallThrough();
    spyOn(svg, 'select').andCallThrough();
    spyOn(svg, 'append').andCallThrough();
    spyOn(svg, 'attr').andCallThrough();
    spyOn(svg, 'text').andCallThrough();
    spyOn(svg, 'data').andCallThrough();
    spyOn(svg, 'enter').andCallThrough();
    spyOn(svg, 'remove').andCallThrough();
    spyOn(svg, 'classed').andCallThrough();
    spyOn(svg, 'on').andCallThrough();
    spyOn(svg, 'call').andCallThrough();
    return svg;
  }

  it('create', function() {
    var chart = Ext.create('Esc.d3.MSpectra', {
      width: 500,
      height: 400
    });
    expect(chart.cutoffCls).toEqual('cutoffline');
    expect(chart.selectedPeakCls).toEqual('selected');
    expect(chart.selectedPeak).toBeFalsy();
    expect(chart.cutoff).toEqual(0);
  });

  it('initScales', function() {
    var chart = Ext.create('Esc.d3.MSpectra', {
      width: 500,
      height: 400,
      data: data,
      axesPadding: [0, 0, 0, 0]
    });
    // mock initSvg
    spyOn(chart, 'getWidth').andReturn(500);
    spyOn(chart, 'getHeight').andReturn(400);

    chart.initScales();

    expect(chart.chartWidth, 500);
    expect(chart.chartHeight, 400);
    expect(chart.ranges).toEqual({
      x: {
        min: 0,
        max: 3
      },
      y: {
        min: 0,
        max: 4
      }
    });
    expect(chart.scales.x.domain()).toEqual([0, 3]);
    expect(chart.scales.x.range()).toEqual([0, 500]);
    expect(chart.scales.y.domain()).toEqual([0, 4]);
    expect(chart.scales.y.range()).toEqual([400, 0]);
  });

  it('initAxes', function() {
    var chart = Ext.create('Esc.d3.MSpectra', {
      width: 500,
      height: 400,
      data: data,
      axesPadding: [0, 0, 0, 0]
    });
    // mock initSvg
    spyOn(chart, 'getWidth').andReturn(500);
    spyOn(chart, 'getHeight').andReturn(400);
    chart.initScales();

    chart.initAxes();

    expect(chart.axes.x.scale()).toEqual(chart.scales.x);
    expect(chart.axes.x.ticks()).toEqual({
      0: chart.ticks.x
    });
    expect(chart.axes.x.orient()).toEqual('bottom');
    expect(chart.axes.y.scale()).toEqual(chart.scales.y);
    expect(chart.axes.y.ticks()).toEqual({
      0: chart.ticks.y
    });
    expect(chart.axes.y.orient()).toEqual('left');
  });

  describe('draw', function() {
    it('!markers + !cutoff', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        data: data,
        axesPadding: [0, 0, 0, 0]
      });
      // mock initSvg
      spyOn(chart, 'getWidth').andReturn(500);
      spyOn(chart, 'getHeight').andReturn(400);
      chart.svg = mockSvg();
      spyOn(chart, 'onMarkersReady');
      chart.initScales();

      chart.draw();

      expect(chart.onMarkersReady).not.toHaveBeenCalled();
      expect(chart.svg.attr).not.toHaveBeenCalledWith('class', chart.cutoffCls);
      expect(chart.svg.attr).toHaveBeenCalledWith('class', 'x axis');
      expect(chart.svg.text).toHaveBeenCalledWith('M/z');
      expect(chart.svg.attr).toHaveBeenCalledWith('class', 'y axis');
      expect(chart.svg.attr).toHaveBeenCalledWith('class', 'mspeak');
    });

    it('markers + cutoff', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        data: data,
        cutoff: 3,
        markers: [{
          mz: 3
        }],
        axesPadding: [0, 0, 0, 0]
      });
      // mock initSvg
      spyOn(chart, 'getWidth').andReturn(500);
      spyOn(chart, 'getHeight').andReturn(400);
      chart.svg = mockSvg();
      spyOn(chart, 'onMarkersReady');
      chart.initScales();

      chart.draw();

      expect(chart.onMarkersReady).toHaveBeenCalled();
      expect(chart.svg.attr).toHaveBeenCalledWith('class', chart.cutoffCls);
    });
  });

  it('undraw', function() {
    var chart = Ext.create('Esc.d3.MSpectra', {
      width: 500,
      height: 400,
      data: data,
      cutoff: 3,
      markers: [{
        mz: 3
      }],
      axesPadding: [0, 0, 0, 0]
    });
    // mock initSvg
    spyOn(chart, 'getWidth').andReturn(500);
    spyOn(chart, 'getHeight').andReturn(400);
    chart.svg = mockSvg();
    spyOn(chart, 'clearPeakSelection');

    chart.undraw();

    expect(chart.clearPeakSelection).toHaveBeenCalled();
    expect(chart.svg.remove).toHaveBeenCalled();
    expect(chart.svg.remove.callCount).toBeGreaterThan(3);
  });

  describe('onToggleMarker', function() {
    it('selectpeak', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        data: data,
        cutoff: 3,
        markers: [{
          mz: 3
        }],
        listeners: {
          'selectpeak': function(mz) {
            // TODO test if event was fired
            expect(mz).toEqual(3);
          }
        }
      });
      spyOn(chart, 'markerSelect').andReturn(3);
      chart.onToggleMarker(3);
      expect(chart.markerSelect).toHaveBeenCalled();
      expect(chart.selectedpeak).toEqual(3);
    });

    it('unselectpeak', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        data: data,
        cutoff: 3,
        markers: [{
          mz: 3
        }],
        listeners: {
          'unselectpeak': function(mz) {
            // TODO test if event was fired
            expect(mz).toEqual(3);
          }
        }
      });
      spyOn(chart, 'markerSelect').andReturn(3);
      chart.selectedpeak = 3;
      chart.onToggleMarker(3);
      expect(chart.markerSelect).toHaveBeenCalled();
      expect(chart.selectedpeak).toBeFalsy();
    });
  });

  it('selectPeak', function() {
    var chart = Ext.create('Esc.d3.MSpectra', {
      width: 500,
      height: 400,
      data: data,
      cutoff: 3,
      markers: [{
        mz: 3
      }]
    });
    // mock initSvg
    chart.chartWidth = 500;
    chart.chartHeight = 400;
    chart.svg = mockSvg();
    spyOn(chart, 'markerSelect');

    chart.selectPeak(3);

    expect(chart.selectedpeak).toEqual(3);
  });

  it('selectPeak when peak already selected', function() {
    var chart = Ext.create('Esc.d3.MSpectra', {
      width: 500,
      height: 400,
      data: data,
      cutoff: 3,
      markers: [{
        mz: 3
      }, {
        mz: 4
      }]
    });
    // mock initSvg
    chart.chartWidth = 500;
    chart.chartHeight = 400;
    chart.svg = mockSvg();
    spyOn(chart, 'markerSelect');

    chart.selectPeak(3);

    expect(chart.selectedpeak).toEqual(3);
  });


  it('clearPeakSelection', function() {
    var chart = Ext.create('Esc.d3.MSpectra', {
      width: 500,
      height: 400,
      data: data,
      cutoff: 3,
      markers: [{
        mz: 3
      }]
    });
    // mock initSvg
    chart.chartWidth = 500;
    chart.chartHeight = 400;
    chart.svg = mockSvg();
    spyOn(chart, 'markerSelect');

    chart.clearPeakSelection();

    expect(chart.selectedpeak).toBeFalsy();
    expect(chart.markerSelect).toHaveBeenCalledWith(false);
  });

  it('setMarkers', function() {
    var markers = [{
      mz: 3
    }];
    var chart = Ext.create('Esc.d3.MSpectra', {
      width: 500,
      height: 400,
      data: data,
      cutoff: 3
    });
    // mock initSvg
    chart.chartWidth = 500;
    chart.chartHeight = 400;
    chart.svg = mockSvg();
    spyOn(chart, 'clearPeakSelection');
    spyOn(chart, 'onMarkersReady');

    chart.setMarkers(markers);

    expect(chart.clearPeakSelection).toHaveBeenCalled();
    expect(chart.svg.selectAll).toHaveBeenCalledWith('.marker');
    expect(chart.svg.remove).toHaveBeenCalled();
    expect(chart.markers).toEqual(markers);
    expect(chart.onMarkersReady).toHaveBeenCalled();
  });

  describe('hasMarkers', function() {
    it('true', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        data: data,
        cutoff: 3,
        markers: [{
          mz: 3
        }]
      });
      expect(chart.hasMarkers()).toBeTruthy();
    });

    it('false', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        data: data,
        cutoff: 3
      });
      expect(chart.hasMarkers()).toBeFalsy();
    });
  });

  describe('onMarkersReady', function() {
    it('nodata', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        cutoff: 3,
        markers: [{
          mz: 3
        }]
      });
      // mock initSvg
      chart.chartWidth = 500;
      chart.chartHeight = 400;
      chart.svg = mockSvg();

      chart.onMarkersReady();

      expect(chart.svg.attr).not.toHaveBeenCalled();
    });

    it('withdata', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        data: data,
        cutoff: 3,
        markers: [{
          mz: 3
        }]
      });
      // mock initSvg
      chart.chartWidth = 500;
      chart.chartHeight = 400;
      chart.svg = mockSvg();

      chart.onMarkersReady();

      expect(chart.svg.attr).toHaveBeenCalledWith('class', 'marker lowermarker annotated');
      expect(chart.svg.attr).toHaveBeenCalledWith('class', 'marker uppermarker annotated');
      expect(chart.svg.attr).toHaveBeenCalledWith('transform', jasmine.any(Function));
      expect(chart.svg.text).toHaveBeenCalledWith(jasmine.any(Function));
      expect(chart.svg.on).toHaveBeenCalledWith('click', jasmine.any(Function));
    });
  });

  describe('onZoom', function() {
    it('no markers', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        data: data,
        axesPadding: [0, 0, 0, 0]
      });
      // mock initSvg
      spyOn(chart, 'getWidth').andReturn(500);
      spyOn(chart, 'getHeight').andReturn(400);
      chart.svg = mockSvg();
      chart.initScales();
      chart.initAxes();
      chart.initZoom();

      chart.onZoom();

      expect(chart.svg.attr).not.toHaveBeenCalledWith('transform', jasmine.any(Function));
      expect(chart.svg.attr).toHaveBeenCalledWith('x1', jasmine.any(Function));
      expect(chart.svg.attr).toHaveBeenCalledWith('x2', jasmine.any(Function));
    });

    it('with markers', function() {
      var chart = Ext.create('Esc.d3.MSpectra', {
        width: 500,
        height: 400,
        data: data,
        axesPadding: [0, 0, 0, 0],
        markers: [{
          mz: 3
        }]
      });
      // mock initSvg
      spyOn(chart, 'getWidth').andReturn(500);
      spyOn(chart, 'getHeight').andReturn(400);
      chart.svg = mockSvg();
      chart.initScales();
      chart.initAxes();
      chart.initZoom();

      chart.onZoom();

      expect(chart.svg.attr).toHaveBeenCalledWith('transform', jasmine.any(Function));
    });
  });
});
