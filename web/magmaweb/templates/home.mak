<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Start</title>
<link rel="stylesheet"
	href="${request.extjsroot}/resources/css/ext-all.css" type="text/css"></link>
<link rel="stylesheet"
	href="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.9.2.custom.css')}"
	type="text/css"></link>
<link rel="stylesheet"  href="${request.extjsroot}/examples/ux/css/ItemSelector.css" type="text/css"></link>
<link rel="stylesheet"  href="${request.extjsroot}/examples/writer/writer.css" type="text/css"></link>
<link rel="stylesheet" href="${request.static_url('magmaweb:static/style.css')}" type="text/css"/>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb-libs.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.9.2.custom.min.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/ChemDoodleWeb-sketcher.js')}"></script>
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
  'Ext.form.Panel',
  'Ext.container.ButtonGroup',
  'Ext.toolbar.Spacer',
  'Ext.container.Viewport',
  'Ext.layout.container.Border',
  'Esc.magmaweb.controller.Scans',
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
    trackResetOnLoad: true,
    layout: {
    	type: 'column'
    },
    defaults: {
    	columnWidth: 0.5,
    	bodyPadding: 5,
    	border: false
    },
    items:[{
    	items: [{
        	title: 'Molecules',
        	margin: '0 0 10 0',
            xtype : 'addstructurefieldset'
        }, {
        	title: 'MS Data',
            xtype : 'uploadmsdatafieldset'
    	}]
    }, {
    	items: [{
    	    xtype : 'metabolizefieldset',
 	        margin: '0 0 10 0'
 	    }, {
    	    xtype : 'annotatefieldset'
    	}]
    }],
    buttons: [{
      text: 'Start from scratch',
      tooltip: 'Do no fill form, but go straight to empty result page which can be filled later',
      href: '${request.route_url('jobfromscratch')}',
      hrefTarget: '_self'
    },{
      text: 'Submit',
      handler: function(){
          if(form.isValid()){
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
                  }
              });
          }
      }
    }, {
      text: 'Reset',
      handler: function() {
          this.up('form').getForm().reset();
          this.up('form').getForm().load({
              url: '${request.route_url('defaults.json')}',
              method: 'GET',
              waitMsg: 'Fetching default settings'
          });
      }
    }]
  });

  form.load({
      url: '${request.route_url('defaults.json')}',
      method: 'GET',
      waitMsg: 'Fetching defaults',
      failure: function(form, action) {
          Ext.Error.raise(action.response.responseText);
      }
  });
  // change settings when tree ms data format is chosen.
  var ms_data_format_combo = form.down('component[name=ms_data_format]');
  scan_controller = Ext.create('Esc.magmaweb.controller.Scans');
  scan_controller.getUploadForm = function() {
    return form;
  };
  scan_controller.application = {};
  scan_controller.application.runInfoUrl = function() {
      return '${request.route_path('defaults.json')}';
  };
  ms_data_format_combo.addListener('change', scan_controller.changeMsDataFormat);
  // hook up example action
  form.down('component[action=loadmsdataexample]').addListener('click', scan_controller.loadExample, scan_controller);
  form.down('component[action=loadmsdataexample2]').addListener('click', scan_controller.loadExample2, scan_controller);

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
        disabled: true
      }, {
        text: 'Help',
        href: "${request.route_url('help')}"
      }, {
  	    text: 'Upload result',
        tooltip: 'Upload a result db for viewing',
        href: "${request.route_url('uploaddb')}",
        hrefTarget: '_self'
	  }, {
        text: 'Workspace',
        tooltip: 'My settings and jobs',
        href: "${request.route_url('workspace')}",
        hrefTarget: '_self'
      }, {
        text: 'Logout',
        hidden: true,
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
	  contentEl: 'logos'
	}]
  };

  Ext.create('Ext.container.Viewport', {
    layout: 'border',
    items: [header, form]
  });

  /**
   * Don't allow metabolization when Molecules tab 'Database' is selected.
   */
  function jobtype_toggler() {
      var molecule_tab_title = form.down('addstructurefieldset').down('tabpanel').getActiveTab().title;
      var ms_format = form.getForm().findField('ms_data_format').getValue();
      var scan = form.getForm().findField('scan');
      var new_visible = false;

      if (molecule_tab_title !== 'Database') {
          new_visible = true;
      }

      var metabolize_fields = form.down('metabolizefieldset');
      var current_visible = metabolize_fields.isVisible();

      if (current_visible !== new_visible) {
          if (new_visible) {
              metabolize_fields.show();
          } else {
              // when toggling from Upload/Draw back to Database the metabolization should be disabled to prevent metabolization with db selected
              metabolize_fields.down('component[name=metabolize]').setValue(false);
              metabolize_fields.hide();
          }
      }

      if (ms_format == 'mzxml' && molecule_tab_title === 'Database') {
          // molecules from database: only one spectral tree allowed
          scan.allowBlank = false;
      } else {
          // molecules from upload: allow full LC-MSn
          // or tree format + database -> no scan needed
          scan.allowBlank = true;
      }
  }

  function applyRole(features) {
      if (features.authenticated || features.anonymous) {
          Ext.ComponentQuery.query('component[text=Login]')[0].hide();
          if (!features.anonymous) {
            // non-anonymous authenticated can logout
            Ext.ComponentQuery.query('component[text=Logout]')[0].show();
          }
      } else {
          Ext.ComponentQuery.query('component[text=Workspace]')[0].hide();
      }
      if (!features.run) {
          // hide interactive job run starting points
          Ext.ComponentQuery.query('component[text=Start from scratch]')[0].hide();
          Ext.ComponentQuery.query('component[text=Upload result]')[0].hide();
      }
      if (features.restricted) {
          Ext.ComponentQuery.query('checkbox[name=metabolize]')[0].boxLabelEl.setHTML('(Only first molecule will be metabolized)');
          form.down('addstructurefieldset').down('tabpanel').addListener('tabchange', jobtype_toggler);
          form.getForm().findField('ms_data_format').addListener('change', jobtype_toggler);
          jobtype_toggler();
      }
      scan_controller.application.features = features;
      scan_controller.applyRole();
  }
  <%!
  import json
  %>
  var features = {
      // should run buttons (upload ms data, upload molecules, metabolize (one molecule), annotate) be shown
      run: false,
      // should logout button be shown
      authenticated: ${json.dumps(request.user is not None)},
      // should login button be hidden and workspace be shown
      anonymous: ${json.dumps(request.registry.settings.get('auto_register', False))|n},
      // should restrictions be applied ie force one spectral tree and metabolize toggle
      restricted: ${json.dumps(request.registry.settings.get('restricted', False))|n}
  };
  applyRole(features);

});
</script>
</head>
<body>
	<div id="sketcher_content" class="x-hidden">
		<script type="text/javascript">
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
	<%include file="logos.mak"/>
</body>
</html>
