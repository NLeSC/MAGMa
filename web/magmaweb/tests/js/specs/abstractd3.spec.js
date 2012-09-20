describe('Esc.d3.Abstract', function() {
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
});