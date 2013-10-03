<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Job status</title>
<link rel="stylesheet"
    href="${request.extjsroot}/resources/css/ext-all.css" type="text/css"></link>
<link rel="stylesheet" href="${request.static_url('magmaweb:static/style.css')}" type="text/css"/>
<style type="text/css">
.status {
  font-size: 200%;
}
</style>
<script type="text/javascript" src="${request.extjsroot}/ext.js"></script>
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
## Comment out below for development or when running sencha build, every class is loaded when commented out
<script type="text/javascript" src="${request.static_url('magmaweb:static/app/resultsApp-all.js')}"></script>
<script type="text/javascript">

Ext.require([
    'Ext.Template',
    'Ext.container.Viewport',
    'Ext.layout.container.Border',
    'Ext.toolbar.Spacer',
    'Ext.container.ButtonGroup',
    'Ext.button.Button',
    'Ext.window.MessageBox',
    'Ext.ProgressBar',
    'Ext.JSON',
    'Ext.Ajax'
]);

Ext.onReady(function() {
  Ext.QuickTips.init();

  var interval = 2000;
  var template = new Ext.Template('<div>Job status:</div><div class="status">{0}</div>');
  template.compile();

  var progress = Ext.create('Ext.ProgressBar', {
      width: '100%',
      border: true,
      listeners: {
          'update': function() {
            if (!progress.isWaiting()) return;
            Ext.Ajax.request({
              url: '${request.route_url('status.json',jobid=jobid)}',
              method: 'GET',
              success: function(response) {
                var result = Ext.JSON.decode(response.responseText);
                if (result.is_complete) {
                    window.location = '${request.route_url('results',jobid=jobid)}';
                } else {
                    Ext.getCmp('status_cont').update([result.status]);
                }
              },
              failure: function(response) {
                  progress.reset();
                  Ext.Msg.alert('Failed to fetch job status', 'Reload page to try again');
              }
            });
          }
      }
  });
  progress.wait({interval: interval});

  var cancel = Ext.create('Ext.button.Button', {
    text: 'Cancel',
    handler: function() {
      // job could complete and redirect while cancelling so stop polling
      progress.reset();
      Ext.Msg.confirm(
        'Cancel job',
        'Are you sure you want cancel job?',
        function(button) {
          if (button === 'yes') {
              Ext.Ajax.request({
                url: '${request.route_url('results',jobid=jobid)}',
                method: 'DELETE',
                success: function() {
                  window.location = '${request.route_url('home')}';
                },
                failure: function() {
                  Ext.Msg.alert('Failed to cancel job', 'Failed to cancel job');
                }
              });
          } else {
             // continue polling when job cancel was cancelled
             progress.wait({interval: interval});
          }
        }
      );
    }
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
        href: "${request.route_url('workspace')}",
        hrefTarget: '_self'
      }, {
        text: 'Logout',
        href: "${request.route_url('logout')}",
        hrefTarget: '_self'
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
      layout: {
          type: 'vbox',
          pack: 'center'
      },
      items: [{
          xtype: 'component',
          tpl: template,
          data: ["${status|n}"],
          id: 'status_cont',
        },
        progress,
        cancel
      ],
      border: false,
      bodyPadding: 50,
      autoScroll: true,
      defaults: {margin: 5}
    }]
  });

  <%!
  import json
  %>
  var authenticated = ${json.dumps(request.user is not None)};
  var anonymous = ${json.dumps(request.registry.settings.get('auto_register', False))|n};
  function user_authenticated(authenticated, anonymous) {
      if (authenticated && anonymous) {
          // anonymous authenticated can not logout
          Ext.ComponentQuery.query('component[text=Logout]')[0].hide();
      }
  };

  user_authenticated(authenticated, anonymous);
});
</script>
</head>
<body>
<%include file="logos.mak"/>
</body>
</html>
