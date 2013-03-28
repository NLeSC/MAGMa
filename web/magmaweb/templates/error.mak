<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Error</title>
<link rel="stylesheet"
  href="${request.extjsroot}/resources/css/ext-all.css" type="text/css"></link>
<link rel="stylesheet" href="${request.static_url('magmaweb:static/style.css')}" type="text/css"/>
<script type="text/javascript" src="${request.extjsroot}/ext.js"></script>
<style type="text/css">

.icon-delete {
  background-image:
    url(${request.extjsroot}/examples/writer/images/delete.png) !important;
}
</style>
<script type="text/javascript">
Ext.Loader.setConfig({
  enabled: true,
//  disableCaching: false, // uncomment to use firebug breakpoints
  paths: {
    'Esc.magmaweb': '${request.static_url('magmaweb:static/app')}',
    'Esc': '${request.static_url('magmaweb:static/esc')}',
    'Ext.ux': '${request.extjsroot}/examples/ux'
  }
});

</script>
## Comment out below for development or when running sencha build, every
## class is loaded when commented out
<script type="text/javascript"
  src="${request.static_url('magmaweb:static/app/resultsApp-all.js')}"></script>
<script type="text/javascript">

Ext.require([
  'Ext.container.Viewport',
  'Ext.layout.container.Border',
  'Ext.toolbar.Spacer',
  'Ext.container.ButtonGroup',
  'Ext.form.Panel',
  'Ext.grid.Panel',
  'Ext.grid.column.Date',
  'Ext.grid.column.Boolean',
  'Ext.grid.column.Action',
  'Ext.grid.plugin.CellEditing',
  'Ext.data.proxy.Rest',
  'Ext.window.MessageBox'
]);

Ext.onReady(function() {
  Ext.QuickTips.init();

  var infoWindow = Ext.create('Ext.window.Window', {
      title: 'Information',
      width: 600,
      autoHeight: true,
      closeAction: 'hide',
      contentEl: 'resultsinfo',
      tools: [{
          type: 'save',
          tooltip: 'Save log file',
          handler: function() {
              window.open('${request.route_url('stderr.txt', jobid=exception.job.id)}', 'Log');
          }
      }]
  });

  var header = {
    border: false,
    region: 'north',
    layout: {
      type: 'hbox',
      align: 'middle',
      padding: 2
    },
    items: [{
      xtype: 'buttongroup',
      columns: 2,
      items: [{
        text: 'Home',
        href: "${request.route_url('home')}",
        hrefTarget: '_self'
      }, {
        text: 'Help',
        href: "${request.route_url('help')}"
      }, {
        text: 'Workspace',
        tooltip: 'My settings and jobs',
        href: "${request.route_url('workspace')}"
      }, {
        text: 'Logout',
        href: "${request.route_url('logout')}",
        hrefTarget: '_self'
      }, {
        text: 'Information',
        tooltip: 'Information about input parameters',
        handler: function() {
            infoWindow.show();
        }
      }]
    }, {
      xtype:'tbspacer',
      flex:1 // aligns buttongroup right
    }, {
      xtype: 'component',
      cls: 'x-title',
      html: '<a href="." data-qtip="<b>M</b>s <b>A</b>nnotation based on in silico <b>G</b>enerated <b>M</b>et<b>a</b>bolites">MAGMa</a>'
    }, {
      xtype:'tbspacer',
      flex:1 // aligns buttongroup right
    }, {
      xtype: 'component',
      contentEl: 'logos'
    }]
  };

  Ext.create('Ext.container.Viewport', {
    layout: 'border',
    items: [header, {
      region: 'center',
      items: [{
        title: 'Calculation failed',
        contentEl: 'error'
      }],
      border: false,
      bodyPadding: 5,
      autoScroll: true,
      defaults: { bodyPadding: 5 }
    }]
  });

  function user_authenticated(authenticated, anonymous) {
      if (authenticated && anonymous) {
          // anonymous authenticated can not logout
          Ext.ComponentQuery.query('component[text=Logout]')[0].hide();
      }
    };

    <%!
    import json
    %>
    var authenticated = ${json.dumps(request.user is not None)};
    var anonymous = ${json.dumps(request.registry.settings.get('auto_register', False))|n};
    user_authenticated(authenticated, anonymous);
});
</script>
</head>
<body>
<%include file="runinfo.mak"/>
<%include file="logos.mak"/>
<div id="error" class="x-hidden">
<h1>
${exception.message}
</h1>
<div>
Please <a href="${request.route_path('home')}">try again</a> with different parameters
</div>
</div>
</body>
</html>
