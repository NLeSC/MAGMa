describe('Esc.d3.Abstract', function() {

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

  var data = [{x:1,y:2},{x:3,y:4}];
  describe('create', function() {
    it('default', function() {
      var chart = Ext.create('Esc.d3.Abstract', {
        width: 400, height: 500
      });
      expect(chart.axesPadding).toEqual([10, 14, 38, 80]);
      expect(chart.ticks).toEqual({x:10, y:4});
      expect(chart.data).toEqual([]);
      expect(chart.hasData()).toBeFalsy();
    });

    it('with data', function() {
      var chart = Ext.create('Esc.d3.Abstract', {
        width: 400, height: 500,
        data: data
      });
      expect(chart.data).toEqual(data);
      expect(chart.hasData()).toBeTruthy();
    });
  });

  describe('onRender', function() {
    it('default', function() {
      var chart = Ext.create('Esc.d3.Abstract', {
        width: 400, height: 500
      });
      spyOn(chart, 'callParent');
      spyOn(chart, 'initSvg')
      chart.body = {
          getWidth: function() { return 400 },
          getHeight: function() { return 500 },
      };

      chart.onRender();

      expect(chart.callParent).toHaveBeenCalled();
      expect(chart.initSvg).toHaveBeenCalled();
    });
  });

  describe('setData', function() {
    it('filled', function() {
      var chart = Ext.create('Esc.d3.Abstract', {
        width: 400, height: 500
      });
      spyOn(chart, 'onDataReady');
      spyOn(chart, 'onDataEmpty');
      chart.setData(data);
      expect(chart.data).toEqual(data);
      expect(chart.hasData()).toBeTruthy();
      expect(chart.onDataReady).toHaveBeenCalled();
      expect(chart.onDataEmpty).not.toHaveBeenCalled();
    });

    it('empty', function() {
      var chart = Ext.create('Esc.d3.Abstract', {
        width: 400, height: 500
      });
      spyOn(chart, 'onDataReady');
      spyOn(chart, 'onDataEmpty');
      chart.setData([]);
      expect(chart.data).toEqual([]);
      expect(chart.hasData()).toBeFalsy();
      expect(chart.onDataReady).not.toHaveBeenCalled();
      expect(chart.onDataEmpty).toHaveBeenCalled();
    });
  });

  describe('onDataEmpty', function() {
    it('!emptyText', function() {
      var chart = Ext.create('Esc.d3.Abstract', {
        width: 400, height: 500
      });
      chart.svg = { append: function() {}};
      spyOn(chart.svg,'append');
      spyOn(chart, 'undraw');
      chart.onDataEmpty();

      expect(chart.undraw).toHaveBeenCalled();
      expect(chart.svg.append).not.toHaveBeenCalled();
    });

    it('emptyText', function() {
      var emptyText = 'Select fragment';
      var chart = Ext.create('Esc.d3.Abstract', {
        width: 400, height: 500, emptyText: emptyText
      });
      chart.svg = {
        append: function() { return this; },
        attr: function() { return this; },
        text: function() { return this; }
      };
      spyOn(chart.svg, 'append').andCallThrough();
      spyOn(chart.svg, 'attr').andCallThrough();
      spyOn(chart.svg, 'text').andCallThrough();
      spyOn(chart, 'undraw');
      chart.chartWidth = 300;
      chart.chartHeight = 400;
      chart.onDataEmpty();
      expect(chart.undraw).toHaveBeenCalled();
      expect(chart.svg.append).toHaveBeenCalledWith('svg:text');
      expect(chart.svg.text).toHaveBeenCalledWith(emptyText);
    });
  });

  it('onDataReady', function() {
    var chart = Ext.create('Esc.d3.Abstract', {
      width: 400, height: 500
    });
    spyOn(chart, 'initScales');
    spyOn(chart, 'undraw');
    spyOn(chart, 'draw');
    chart.onDataReady();
    expect(chart.initScales).toHaveBeenCalled();
    expect(chart.undraw).toHaveBeenCalled();
    expect(chart.draw).toHaveBeenCalled();
  });

  it('initScales', function() {
      var chart = Ext.create('Esc.d3.Abstract', {
          width: 400, height: 500,
          axesPadding: [20, 10, 40, 80]
      });
      // mock initSvg
      spyOn(chart, 'getWidth').andReturn(400);
      spyOn(chart, 'getHeight').andReturn(500);

      chart.initScales();

      expect(chart.chartWidth, 400 - 10 - 40);
      expect(chart.chartHeight, 500 - 20 - 80);
  });

  it('resetScales', function() {
    var chart = Ext.create('Esc.d3.Abstract', {
      width: 400, height: 500
    });
    chart.ranges = { x: {min: 0, max: 100}, y: {min: 200, max: 300} };
    chart.scales.x = { domain: function() { return this; } };
    chart.scales.y = { domain: function() { return this; } };
    spyOn(chart.scales.x, 'domain');
    spyOn(chart.scales.y, 'domain');
    spyOn(chart, 'initZoom');
    spyOn(chart, 'onZoom');

    chart.resetScales();

    expect(chart.scales.x.domain).toHaveBeenCalledWith([0, 100]);
    expect(chart.scales.y.domain).toHaveBeenCalledWith([200, 300]);
    expect(chart.initZoom).toHaveBeenCalled();
    expect(chart.onZoom).toHaveBeenCalled();
  });

  it('zoomBehavior', function() {
      var chart = Ext.create('Esc.d3.Abstract', {
      width: 400, height: 500,
      scales: {
        x: d3.scale.linear(),
        y: d3.scale.linear()
      }
      });

      var zoom = chart.zoomBehavior();

      expect(chart.scales0.x).toBeDefined();
      expect(chart.scales0.y).toBeDefined();
  });

  it('enableZoom', function() {
	  var chart = Ext.create('Esc.d3.Abstract', {
	    width: 400, height: 500,
	    zoom: {x: true, y: false}
	  });
	  spyOn(chart, 'initZoom');

	  chart.setZoom('y', true);

	  expect(chart.zoom.y).toEqual(true);
	  expect(chart.initZoom).toHaveBeenCalledWith();
  });

  describe('onZoom', function() {
      var chart;

	  beforeEach(function() {
		chart = Ext.create('Esc.d3.Abstract', {
		  width: 400, height: 500,
		  zoom: {x: true, y: true}
		});
		// mock initSvg
		spyOn(chart, 'getWidth').andReturn(500);
		spyOn(chart, 'getHeight').andReturn(400);
		chart.svg = mockSvg();
		chart.initScales();
		chart.scales.x = d3.scale.linear().domain([1,3]).range([0, chart.chartWidth]);
		chart.scales.y = d3.scale.linear().domain([2,4]).range([chart.chartHeight, 0]);
		chart.initZoom();
	  });

	  it('Zoom x,y without d3.event', function() {
	      // perform a zoom/pan
	      chart.scales.x.domain([5,6]);
	      chart.scales.y.domain([7,8]);

	      chart.onZoom();

	      // expect scales to reset back
	      expect(chart.scales.x.domain()).toEqual([1,3]);
	      expect(chart.scales.y.domain()).toEqual([2,4]);
	      // expect svg axis to draw
	      expect(chart.svg.select).toHaveBeenCalledWith('.x.axis');
	      expect(chart.svg.select).toHaveBeenCalledWith('.y.axis');
	  });

	  it('Zoom x with d3.event', function() {
		  chart.zoom.y = false;
		  d3.event = {
		      translate: [9,10],
		      scale: 2
		  };

	      chart.onZoom();

	      expect(chart.scales.x.domain()).toEqual([0.9778325123152709, 1.977832512315271]);
	      expect(chart.scales.y.domain()).toEqual([2,4]);
	      expect(chart.svg.select).toHaveBeenCalledWith('.x.axis');
	      expect(chart.svg.select).not.toHaveBeenCalledWith('.y.axis');
	  });

	  it('Zoom y with d3.event', function() {
		  chart.zoom.x = false;
		  d3.event = {
		      translate: [9,10],
		      scale: 2
		  };

	      chart.onZoom();

	      expect(chart.scales.x.domain()).toEqual([1,3]);
	      expect(chart.scales.y.domain()).toEqual([2,3]);
	      expect(chart.svg.select).not.toHaveBeenCalledWith('.x.axis');
	      expect(chart.svg.select).toHaveBeenCalledWith('.y.axis');
	  });
  });
});
