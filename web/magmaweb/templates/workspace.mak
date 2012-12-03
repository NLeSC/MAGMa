<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
<title>MAGMa - Workspace</title>
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

  var form = Ext.create('Ext.form.Panel', {
    title: 'User',
    autoScroll: true,
    items:[{
    	xtype: 'displayfield',
    	fieldLabel: 'User id',
    	value: '${request.user.userid}'
    }, {
        xtype: 'textfield',
    	fieldLabel: 'Name',
    	width: 400,
    	value: '${request.user.displayname}'
    }, {
        xtype: 'textfield',
    	fieldLabel: 'Email',
    	width: 400,
    	value: '${request.user.email}',
    	vtype: 'email'
    }],
    buttons: [{
      text: 'Update'
    }]
  });

  <%!
  import json
  %>
  var job_store = Ext.create('Ext.data.Store', {
     fields: ['id', 'description', 'ms_filename', 'created_at', 'url'],
     data: ${json.dumps(jobs)|n}
  });

  var job_grid = Ext.create('Ext.grid.Panel', {
    title: 'Jobs',
    store: job_store,
    columns: [{
      text: 'ID', dataIndex: 'id', renderer: function(v, m, r) {
        return Ext.String.format('<a href="{0}">{1}</a>', r.data.url, v);
      },
      width: 220
    }, {
      text: 'Description', dataIndex: 'description',
      width: 600
    }, {
      text: 'ms_filename', dataIndex: 'ms_filename'
    }, {
      text: 'created_at', dataIndex: 'created_at'
    }],
    listeners: {
      select: function(rm, r) {
      	window.location = r.data.url;
      }
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
        tooltip: 'Goto help pages',
        disabled: true
      }, {
        text: 'Workspace',
        tooltip: 'My settings and jobs',
        disabled: true
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
      cls: 'x-logo',
      html: '<a href="http://www.esciencecenter.nl"></a>'
    }]
  };

  Ext.create('Ext.container.Viewport', {
    layout: 'border',
    items: [header, {
    	region: 'center',
    	items: [form, job_grid],
        border: false,
    	bodyPadding: 5,
	    defaults: { bodyPadding: 5 },
    }]
  });
});
</script>
</head>
<body>
<div id="welcome" class="x-hidden">
<h1>Workspace</h1>
    <h2>User settings</h2>
    <form>
    <div>
        <label for="userid">User id</label>
        <input name="userid" disabled value="${request.user.userid}"></input>
        </div><div>
        <label for="displayname">Display Name</label>
        <input name="displayname" value="${request.user.displayname}" size=50></input>
        </div><div>
        <label for="email">Email</label>
        <input name="email" value="${request.user.email}" size=50></input>
        </div>
        <button>Update</button>
    </form>
    <h2>Jobs</h2>
    <table>
<thead><tr><th>ID</th><th>Description</th><th>MS filename</th><th>Created at</th></tr></thead>
<tbody>
        % for job in jobs:
        <tr><td><a href="${request.route_url('results',jobid=job['id'])}">${job['id']}</a></td><td><a href="${request.route_url('results',jobid=job['id'])}">${job['description']}</a></td>
        <td>${job['ms_filename']}</td><td>${job['created_at']}</td>
        </tr>
        % endfor
    </tbody></table>

    <hr></hr>
    <a href="${request.route_url('home')}">Home</a>
    </div>
</body>
</html>
