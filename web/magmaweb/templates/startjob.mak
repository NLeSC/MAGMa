<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Ms Annotation based on in silico Generated Metabolites</title>
<link rel="stylesheet"
	href="${request.extjsroot}/resources/css/ext-all.css" type="text/css"></link>
<link rel="stylesheet"
	href="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.9.2.custom.css')}"
	type="text/css"></link>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb-libs.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.9.2.custom.min.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/ChemDoodleWeb-sketcher.js')}"></script>
<script type="text/javascript" src="${request.extjsroot}/ext.js"></script>
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

</script>
## Comment out below for development or when running sencha build, every class is loaded when commented out
<script type="text/javascript" src="${request.static_url('magmaweb:static/app/resultsApp-all.js')}"></script>
<script type="text/javascript">

Ext.require([
  'Ext.form.Panel',
  'Ext.container.ButtonGroup',
  'Ext.toolbar.Spacer',
  'Ext.container.Viewport',
  'Ext.layout.container.Border',
  'Esc.magmaweb.view.scan.UploadFieldSet',
  'Esc.magmaweb.view.metabolite.AddFieldSet',
  'Esc.magmaweb.view.metabolite.MetabolizeFieldSet',
  'Esc.magmaweb.view.fragment.AnnotateFieldSet'
]);

Ext.onReady(function() {
  Ext.QuickTips.init();

  var form = Ext.create('Ext.form.Panel', {
    border: false,
    region: 'center',
    bodyPadding: 5,
    defaults: { bodyPadding: 5 },
    autoScroll: true,
    items:[{
        contentEl: 'welcome',
        border: false
    }, {
    	title: 'Molecules',
        xtype : 'addstructurefieldset'
    }, {
    	title: 'MS Data',
        xtype : 'uploadmsdatafieldset'
    }, {
        xtype : 'metabolizefieldset',
        checkboxToggle: true,
        checkboxName: 'metabolize',
        collapsed : true,
        collapsible : true
    }, {
        xtype : 'annotatefieldset',
        collapsed : true,
        collapsible : true
    }],
    buttons: [{
      text: 'Start from scratch',
      tooltip: 'Do no fill form, but go straight to empty result page which can be filled later',
      href: '${request.route_url('jobfromscratch')}',
      hrefTarget: '_self'
    },{
      text: 'Submit',
      handler: function(){
          var form = this.up('form').getForm();
          var mol = sketcher.getMolecule();
          if (mol.bonds.length > 0) {
              form.setValues({
                 structures_format: 'sdf',
                 structures_area: ChemDoodle.writeMOL(mol)
              });
          }
          if(form.isValid()){
              // TODO test if structures textarea or file is filled
              form.submit({
                  url: '${request.route_url('startjob')}',
                  waitMsg: 'Uploading your data...',
                  submitEmptyText: false,
                  success: function(fp, o) {
                      window.location = '${request.application_url}/status/'+o.result.jobid;
                  },
                  failure: function(form, action) {
                	  if ('msg' in action.result) {
	                	  Ext.Msg.show({
	                		  title: 'Something went wrong submitting job',
	                		  msg: action.result.msg,
	                		  icon: Ext.MessageBox.ERROR,
	                		  buttons: Ext.MessageBox.OK,
	                	  });
                	  }
                      console.log(action.failureType);
                      console.log(action.result);
                  }
              });
          }
      }
    }, {
      text: 'Reset',
      handler: function() {
        this.up('form').getForm().reset();
      }
    }]
  });
  form.load({
      url: '${request.route_url('defaults.json')}',
      method: 'GET',
      waitMsg: 'Fetching defaults'
  });
  // hook up example action
  var example_button = form.down('component[action=loadmsdataexample]');
  example_button.addListener('click', function() {
	  form.load({
	      url: '${request.route_url('defaults.json', _query={'selection': 'example'})}',
	      method: 'GET',
	      waitMsg: 'Fetching example settings'
	  });
  });
  // change settings when tree ms data format is chosen.
  var ms_data_format_combo = form.down('component[name=ms_data_format]');
  ms_data_format_combo.addListener('change', function(field, value) {
	if (value == 'tree') {
		form.getForm().setValues({
		    'ms_intensity_cutoff': 0,
		    'msms_intensity_cutoff': 0,
		    'abs_peak_cutoff': 0
		});
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
	    text: 'Help',
	    tooltip: 'Goto help pages',
	    disabled: true
	  }, {
	    text: 'Upload result',
	    tooltip: 'Upload a result db for viewing',
	    href: '${request.route_url('uploaddb')}',
        hrefTarget: '_self'
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
	<div id="sketcher_content" class="x-hidden">
		<script language="javascript">
			var sketcher = new ChemDoodle.SketcherCanvas(
		        'sketcher_canvas', 500, 300, {
	        		useServices: false, oneMolecule: true
	        	}
		    );
			sketcher.repaint();
			sketcher.toolbarManager.setup();
			sketcher.toolbarManager.buttonSave.disable();
			sketcher.toolbarManager.buttonOpen.disable();
		</script>
	</div>
	<div id="welcome" class="x-hidden">
		<h1>
			Welcome to the <b>M</b>s <b>A</b>nnotation based on in silico <b>G</b>enerated
			<b>M</b>et<b>a</b>bolites application
		</h1>
	</div>
</body>
</html>
