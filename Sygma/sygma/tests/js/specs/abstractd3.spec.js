describe('Ext.esc.AbstractD3', function() {
  var data = [{x:1,y:2},{x:3,y:4}];
  describe('create', function() {
    it('default', function() {
      var chart = Ext.create('Ext.esc.AbstractD3', {
        width: 400, height: 500
      });
      expect(chart.axesPadding).toEqual([16, 5, 38, 80]);
      expect(chart.ticks).toEqual({x:10, y:4});
      expect(chart.data).toEqual([]);
      expect(chart.hasData()).toBeFalsy();
    });

    it('with data', function() {
      var chart = Ext.create('Ext.esc.AbstractD3', {
        width: 400, height: 500,
        data: data
      });
      expect(chart.data).toEqual(data);
      expect(chart.hasData()).toBeTruthy();
    });
  });

  describe('onRender', function() {
    it('default', function() {
      var chart = Ext.create('Ext.esc.AbstractD3', {
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

    it('percentsize', function() {
      var chart = Ext.create('Ext.esc.AbstractD3', {
        width: '50%', height: '30%'
      });
      spyOn(chart, 'callParent');
      spyOn(chart, 'initSvg');
      chart.body = {
          getWidth: function() { return 1 },
          getHeight: function() { return 1 },
      };
      chart.onRender();
      expect(chart.callParent).toHaveBeenCalled();
      expect(chart.initSvg).not.toHaveBeenCalled();
      // mock layout
      chart.body.getWidth = function() { return 500 };
      chart.body.getHeight = function() { return 500 };
      chart.fireEvent('afterlayout', chart);
      expect(chart.initSvg).toHaveBeenCalled();
    });
  });

  it('resize', function() {
    var chart = Ext.create('Ext.esc.AbstractD3', {
      width: 400, height: 500
    });
    spyOn(chart, 'on');
    spyOn(chart, 'initSvg')
    chart.body = {
        getWidth: function() { return 400 },
        getHeight: function() { return 500 },
    };
    var svg = {
        attr: function(key, value) {},
        select: function() { return this; },
    };
    spyOn(svg, 'attr');
    spyOn(d3, 'select').andReturn(svg);
    chart.fireEvent('resize', chart, 600, 700);
    expect(svg.attr).toHaveBeenCalledWith('width', 600);
    expect(svg.attr).toHaveBeenCalledWith('height', 700);
  });

  describe('setData', function() {
    it('filled', function() {
      var chart = Ext.create('Ext.esc.AbstractD3', {
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
      var chart = Ext.create('Ext.esc.AbstractD3', {
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
      var chart = Ext.create('Ext.esc.AbstractD3', {
        width: 400, height: 500
      });
      chart.svg = { append: function() {}};
      spyOn(chart.svg,'append');
      chart.onDataEmpty();
      expect(chart.svg.append).not.toHaveBeenCalled();
    });

    it('emptyText', function() {
      var emptyText = 'Select fragment';
      var chart = Ext.create('Ext.esc.AbstractD3', {
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
      chart.chartWidth = 300;
      chart.chartHeight = 400;
      chart.onDataEmpty();
      expect(chart.svg.append).toHaveBeenCalled();
      expect(chart.svg.attr).toHaveBeenCalledWith('x', 150);
      expect(chart.svg.attr).toHaveBeenCalledWith('y', 200);
      expect(chart.svg.text).toHaveBeenCalledWith(emptyText);
    });
  });

  it('onDataReady', function() {
    var chart = Ext.create('Ext.esc.AbstractD3', {
      width: 400, height: 500
    });
    spyOn(chart, 'initScales');
    spyOn(chart, 'initAxes');
    chart.onDataReady();
    expect(chart.initScales).toHaveBeenCalled();
    expect(chart.initAxes).toHaveBeenCalled();
  });
});