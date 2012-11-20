<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
<title>MAGMa - Ms Annotation based on in silico Generated Metabolites</title>
<link rel="stylesheet"
    href="${request.extjsroot}/resources/css/ext-all.css" type="text/css"></link>
<script type="text/javascript" src="${request.extjsroot}/ext-all.js"></script>
<style type="text/css">
.x-logo a {
  font-size: 40px;
  padding-left: 520px;
  padding-top: 3px; /* aligns app title with text in logo  */
  background: url(${request.static_url('magmaweb:static/ESCIENCE_log_B_nl_long_cyanblack.jpg')}) no-repeat 5px 4px;
}

.x-title a {
  font-size: 40px;
  font-weight: bold;
  color: #333;
  text-decoration:none;
  margin-left: auto;
  margin-right: auto;
  padding-top: 3px; /* aligns app title with text in logo  */
}

#welcome h1 {
    font-size: 200%;
    padding-top: 40px;
    padding-bottom: 40px;
    color: #333;
}
</style>
<script type="text/javascript">
if (!window.console) window.console = {};
if (!window.console.log) window.console.log = function() {};

Ext.Loader.setConfig({
  enabled: true,
//  disableCaching: false, // uncomment to use firebug breakpoints
  paths: {
    'Esc.magmaweb': '${request.static_url('magmaweb:static/app')}',
    'Esc': '${request.static_url('magmaweb:static/esc')}',
    'Ext.ux': '${request.extjsroot}/examples/ux'
  }
});

Ext.require([
  'Esc.magmaweb.view.scan.UploadFieldSet',
  'Esc.magmaweb.view.metabolite.AddFieldSet',
  'Esc.magmaweb.view.metabolite.MetabolizeFieldSet',
  'Esc.magmaweb.view.fragment.AnnotateFieldSet'
]);

Ext.onReady(function() {
  Ext.QuickTips.init();

  var form = Ext.create('Ext.panel.Panel', {
    border: false,
    region: 'center',
    bodyPadding: 5,
    defaults: { bodyPadding: 5 },
    autoScroll: true,
    items:[{
        contentEl: 'welcome',
        border: false
    }],
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
        text: 'Help',
        tooltip: 'Goto help pages',
        disabled: true
      }, {
        text: 'Workspace',
        tooltip: 'My settings and jobs',
        href: "${request.route_url('workspace')}",
        hrefTarget: '_self'
      }, {
        text: 'Logout',
        href: "${request.route_url('logout')}",
        hrefTarget: '_self'
      }, {
        text: 'Login',
        href: "${request.route_url('login')}",
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
      cls: 'x-logo',
      html: '<a href="http://www.esciencecenter.nl"></a>'
    }]
  };

  Ext.create('Ext.container.Viewport', {
    layout: 'border',
    items: [header, form]
  });

  <%!
  import json
  %>
  function user_authenticated(toggle) {
     if (toggle) {
    	 // hide login
    	 Ext.ComponentQuery.query('component[text=Login]')[0].hide();
     } else {
    	 // hide workspace+logout
         Ext.ComponentQuery.query('component[text=Workspace]')[0].hide();
         Ext.ComponentQuery.query('component[text=Logout]')[0].hide();
     }
  };

  user_authenticated(${json.dumps(request.user is not None)});

});
</script>
</head>
<body>
    <div id="welcome" class="x-hidden">
        <h1>
            Welcome to the <b>M</b>s <b>A</b>nnotation based on in silico <b>G</b>enerated
            <b>M</b>et<b>a</b>bolites application
        </h1>
        <ul>
        <li><a href="${request.route_url('startjob')}">Start new calculation</a></li>
        </ul>
    </div>
</body>
</html>
