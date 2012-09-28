describe('Esc.d3.Chromatogram', function() {
  var data = [{rt:1, intensity: 100, id:4}, {rt:2, intensity: 50, id:5}];

  function mockSvg() {
    var svg = {
      selectAll: function() { return this; },
      select: function() { return this; },
      append: function() { return this; },
      attr: function() { return this; },
      call: function() { return this; },
      data: function() { return this; },
      enter: function() { return this; },
      remove: function() { return this; },
      classed: function() { return this; },
      on: function() { return this; },
      text: function() { return this; }
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
    return svg;
  }

  it('create', function() {
    var chart = Ext.create('Esc.d3.Chromatogram', {
      width: 500, height: 400
    });
    expect(chart.cutoffCls).toEqual('cutoffline');
    expect(chart.selectedScanCls).toEqual('selected');
    expect(chart.selectedScan).toEqual(-1);
    expect(chart.cutoff).toEqual(2000000);
  });

  it('initScales', function() {
    var chart = Ext.create('Esc.d3.Chromatogram', {
      width: 500, height: 400, data: data,
      axesPadding: [0, 0, 0, 0]
    });
    // mock initSvg
    spyOn(chart, 'getWidth').andReturn(500);
    spyOn(chart, 'getHeight').andReturn(400);

    chart.initScales();

    expect(chart.ranges).toEqual({
      x: { min:0, max: 2},
      y: { min:0, max: 100}
    });
    expect(chart.scales.x.domain()).toEqual([0,2]);
    expect(chart.scales.x.range()).toEqual([0,500]);
    expect(chart.scales.y.domain()).toEqual([0,100]);
    expect(chart.scales.y.range()).toEqual([400,0]);
    expect(chart.axes.x.scale()).toEqual(chart.scales.x);
    expect(chart.axes.x.ticks()).toEqual({ 0:chart.ticks.x});
    expect(chart.axes.x.orient()).toEqual('bottom');
    expect(chart.axes.y.scale()).toEqual(chart.scales.y);
    expect(chart.axes.y.ticks()).toEqual({ 0:chart.ticks.y});
    expect(chart.axes.y.orient()).toEqual('left');
  });

  describe('draw', function() {
    it('!markers', function() {
      var chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400, data: data, cutoff: 10,
        axesPadding: [0, 0, 0, 0]
      });
      // mock initSvg
      spyOn(chart, 'getWidth').andReturn(500);
      spyOn(chart, 'getHeight').andReturn(400);
      chart.svg = mockSvg();
      spyOn(chart,'onMarkersReady');
      chart.initScales();

      chart.draw();

      expect(chart.onMarkersReady).not.toHaveBeenCalled();
      expect(chart.svg.attr).toHaveBeenCalledWith('class', chart.cutoffCls);
      expect(chart.svg.attr).toHaveBeenCalledWith('class', 'x axis');
      expect(chart.svg.text).toHaveBeenCalledWith('Retention time (min)');
      expect(chart.svg.attr).toHaveBeenCalledWith('class', 'y axis');
      expect(chart.svg.text).toHaveBeenCalledWith('Intensity');
      expect(chart.svg.attr).toHaveBeenCalledWith('class', 'peak');
      expect(chart.svg.attr).toHaveBeenCalledWith('class', 'line');
      expect(chart.svg.data).toHaveBeenCalledWith(data);
    });

    it('markers', function() {
      var chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400, data: data,
        cutoff: 10, markers: [data[0]],
        axesPadding: [0, 0, 0, 0]
      });
      // mock initSvg
      spyOn(chart, 'getWidth').andReturn(500);
      spyOn(chart, 'getHeight').andReturn(400);

      chart.svg = mockSvg();
      spyOn(chart,'onMarkersReady');
      chart.initScales();

      chart.draw();

      expect(chart.onMarkersReady).toHaveBeenCalled();
    });
  });

  it('undraw', function() {
      var chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400, data: data,
        cutoff: 3, markers: [{mz: 3}],
        axesPadding: [0, 0, 0, 0]
      });
      // mock initSvg
      spyOn(chart, 'getWidth').andReturn(500);
      spyOn(chart, 'getHeight').andReturn(400);
      chart.svg = mockSvg();
      spyOn(chart, 'clearScanSelection');

      chart.undraw();

      expect(chart.clearScanSelection).toHaveBeenCalled();
      expect(chart.svg.remove).toHaveBeenCalled();
      expect(chart.svg.remove.callCount).toBeGreaterThan(5);
  });

  describe('onZoom', function() {
      it('no markers', function() {
          var chart = Ext.create('Esc.d3.Chromatogram', {
            width: 500, height: 400, data: data,
            axesPadding: [0, 0, 0, 0]
          });
          // mock initSvg
          spyOn(chart, 'getWidth').andReturn(500);
          spyOn(chart, 'getHeight').andReturn(400);
          chart.svg = mockSvg();
          chart.initScales();

          chart.onZoom();

          expect(chart.svg.select).toHaveBeenCalledWith('.x.axis');
          expect(chart.svg.select).toHaveBeenCalledWith('path.line');
          expect(chart.svg.select).toHaveBeenCalledWith('path.metaboliteline');
          expect(chart.svg.attr).not.toHaveBeenCalledWith('transform', jasmine.any(Function));
          expect(chart.svg.attr).toHaveBeenCalledWith('x1', jasmine.any(Function));
          expect(chart.svg.attr).toHaveBeenCalledWith('x2', jasmine.any(Function));
      });

      it('with markers', function() {
          var chart = Ext.create('Esc.d3.Chromatogram', {
            width: 500, height: 400, data: data,
            axesPadding: [0, 0, 0, 0], markers: [{mz: 3}]
          });
          // mock initSvg
          spyOn(chart, 'getWidth').andReturn(500);
          spyOn(chart, 'getHeight').andReturn(400);
          chart.svg = mockSvg();
          chart.initScales();

          chart.onZoom();

          expect(chart.svg.attr).toHaveBeenCalledWith('transform', jasmine.any(Function));
      });
  });

  describe('selectScan', function() {
    var chart;
    beforeEach(function() {
      chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400, data: data,
        cutoff: 10, markers: data
      });
      // mock initSvg
      chart.chartWidth = 500;
      chart.chartHeight = 400;
      chart.svg = mockSvg();
    });

    it('select1', function() {
      spyOn(chart,'markerSelect');

      chart.selectScan(data[0].id);

      expect(chart.markerSelect).toHaveBeenCalled();
      expect(chart.selectedScan).toEqual(data[0].id);
    });

    it('select0', function() {
      spyOn(chart,'markerSelect');

      chart.selectScan(null);

      expect(chart.markerSelect).toHaveBeenCalled();
      expect(chart.selectedScan).toEqual(-1);
    });

    it('another scan already selected -> deselect', function() {
        chart.selectScan(data[0].id);

        spyOn(chart,'markerSelect');
        spyOn(chart,'fireEvent');

        chart.selectScan(data[1].id);

        expect(chart.fireEvent).toHaveBeenCalledWith('unselectscan', 4);
        expect(chart.selectedScan).toEqual(5);
    });

    it('another scan already selected -> deselect silenced', function() {
        chart.selectScan(data[0].id);

        spyOn(chart,'markerSelect');
        spyOn(chart,'fireEvent');

        chart.selectScan(data[1].id, true);

        expect(chart.fireEvent).not.toHaveBeenCalledWith('unselectscan', 4);
        expect(chart.selectedScan).toEqual(5);
    });

    it('same scan already selected -> no deselect', function() {
        chart.selectScan(data[0].id);

        spyOn(chart,'markerSelect');
        spyOn(chart,'fireEvent');

        chart.selectScan(data[0].id);

        expect(chart.fireEvent).not.toHaveBeenCalledWith('unselectscan', 4);
        expect(chart.selectedScan).toEqual(4);
    });
  });

  it('clearScanSelection', function() {
    var chart = Ext.create('Esc.d3.Chromatogram', {
      width: 500, height: 400, data: data,
      cutoff: 3, markers: data
    });
    // mock initSvg
    chart.chartWidth = 500;
    chart.chartHeight = 400;
    chart.svg = mockSvg();
    spyOn(chart, 'markerSelect');

    chart.clearScanSelection();

    expect(chart.selectedScan).toEqual(-1);
    expect(chart.markerSelect).toHaveBeenCalledWith(false);
  });

  describe('setMarkers', function() {
    it('no scan selected', function() {
      var markers = data;
      var chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400, data: data,
        cutoff: 3
      });
      // mock initSvg
      chart.chartWidth = 500;
      chart.chartHeight = 400;
      chart.svg = mockSvg();
      spyOn(chart, 'clearScanSelection');
      spyOn(chart, 'onMarkersReady');

      chart.setMarkers(markers);

      expect(chart.clearScanSelection).toHaveBeenCalled();
      expect(chart.svg.selectAll).toHaveBeenCalledWith('.marker');
      expect(chart.svg.remove).toHaveBeenCalled();
      expect(chart.markers).toEqual(markers);
      expect(chart.onMarkersReady).toHaveBeenCalled();
    });

    it('scan selected', function() {
      var markers = data;
      var chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400, data: data,
        cutoff: 3, markers: markers, selectedScan: data[0].id
      });
      // mock initSvg
      chart.chartWidth = 500;
      chart.chartHeight = 400;
      chart.svg = mockSvg();
      spyOn(chart, 'clearScanSelection');
      spyOn(chart, 'onMarkersReady');
      spyOn(chart, 'selectScan');

      chart.setMarkers(markers);

      expect(chart.selectScan).toHaveBeenCalledWith(data[0].id);
      expect(chart.selectedScan).toEqual(data[0].id);
    });
  });

  describe('hasMarkers', function() {
    it('true', function() {
      var chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400, data: data,
        cutoff: 3, markers: data
      });
      expect(chart.hasMarkers()).toBeTruthy();
    });

    it('false', function() {
      var chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400, data: data,
        cutoff: 3
      });
      expect(chart.hasMarkers()).toBeFalsy();
    });
  });

  describe('onMarkersReady', function() {
    it('nodata', function() {
      var chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400,
        cutoff: 3, markers: data
      });
      // mock initSvg
      chart.chartWidth = 500;
      chart.chartHeight = 400;
      chart.svg = mockSvg();

      chart.onMarkersReady();

      expect(chart.svg.attr).not.toHaveBeenCalled();
    });

    it('withdata', function() {
      var chart = Ext.create('Esc.d3.Chromatogram', {
        width: 500, height: 400, data: data,
        cutoff: 3, markers: data
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

  it('setExtractedIonChromatogram', function() {
    var chart = Ext.create('Esc.d3.Chromatogram', {
      width: 500, height: 400, data: data,
      axesPadding: [0, 0, 0, 0]
    });
    // mock initSvg
    spyOn(chart, 'getWidth').andReturn(500);
    spyOn(chart, 'getHeight').andReturn(400);

    chart.svg = mockSvg();

    var eic = data;
    chart.initScales(); // line drawing requires line factory which is made there
    chart.setExtractedIonChromatogram(eic);

    expect(chart.metabolitedata).toEqual(eic);
    expect(chart.svg.attr).toHaveBeenCalledWith('class', 'metaboliteline');
  });
});