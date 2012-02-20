<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <title>MAGMa - Homepage</title>
<link rel="stylesheet" href="${request.extjsroot}/resources/css/ext-all.css" type="text/css"></link>
<script type="text/javascript" src="${request.extjsroot}/ext-all.js"></script>
<script type="text/javascript">

Ext.Loader.setConfig({
  enabled: true,
//  disableCaching: false, // uncomment to use firebug breakpoints
  paths: {
    'Esc': '${request.static_url('magmaweb:static/esc')}',
    'Ext.ux': '${request.extjsroot}/examples/ux'
  }
});

Ext.require([
  'Ext.form.Panel',
  'Ext.layout.container.Column',
  'Ext.form.FieldSet',
  'Ext.form.RadioGroup',
  'Ext.form.field.File',
  'Ext.form.field.Number',
  'Ext.form.field.Radio'
]);

Ext.onReady(function() {
  Ext.QuickTips.init();

  var form = Ext.create('Ext.form.Panel', {
    frame: true,
    width: 600,
    title: 'MAGMa submit form',
    bodyPadding: '10 10 0',
    renderTo: Ext.getBody(),
    layout: 'column',
    items:[{
       xtype: 'container',
       items:[{
          xtype: 'textarea',
          fieldLabel: 'Metabolites (one SMILES string per line)',
          name: 'structures',
          height: 300,
        }, {
          fieldLabel: 'MS/MS data (mzxml)',
          name: 'db',
          xtype: 'filefield'
        }]
    }, {
        xtype: 'container',
        items: [{
          xtype: 'fieldset',
          title: 'Generate metabolite options',
          items: [{
              fieldLabel: 'Maximum number of reaction steps',
              name: 'n_reaction_steps',
              xtype: 'numberfield',
              value: 2,
              maxValue: 10,
              minValue: 0,
              decimalPrecision: 0
          },{
              xtype:'combobox',
              fieldLabel: 'Metabolism types',
              store:['phase1', 'phase2'],
              multiSelect: true,
              value: 'phase1, phase2',
              name: 'metabolism_types'
          }]
      },{
          xtype: 'fieldset',
          title: 'MSpectra options',
          items: [{
              xtype: 'numberfield',
              name: 'precursor_mz_precision',
              fieldLabel: 'Precision for matching precursor mz with peak mz in parent scan',
              value: 0.001,
              decimalPrecision: 5
          },{
              xtype: 'numberfield',
              name: 'abs_peak_cutoff',
              fieldLabel: 'Absolute intensity threshold for storing peaks in database',
              value: 1000,
              decimalPrecision: 5
          },{
              xtype: 'numberfield',
              name: 'rel_peak_cutoff',
              fieldLabel: 'Absolute intensity threshold for storing peaks in database',
              value: 0.01,
              decimalPrecision: 5
          }]
      }]
    }, {
        xtype: 'container',
        items: [{
	        xtype: 'fieldset',
	        title: 'Annotate options',
	        items: [{
	            fieldLabel: 'Maximum number of bonds broken in substructures generated from metabolites',
	            name: 'max_broken_bonds',
	            xtype: 'numberfield',
	            value: 4,
	            maxValue: 10,
	            minValue: 0,
	            decimalPrecision: 0
	        },{
	            fieldLabel: 'Ionisation',
	            xtype: 'radiogroup',
	            columns: 2,
	            items: [
	              { boxLabel:'-', checked: true, name:'ionisation', inputValue: -1 },
	              { boxLabel:'+', checked: false, name:'ionisation', inputValue: 1 }
	            ]
	        },{
	            xtype:'checkbox',
	            fieldLabel: 'Use fragmentation',
	            checked: 'checked',
	            name: 'use_fragmentation'
	        },{
	            xtype:'checkbox',
	            fieldLabel: 'Annotate only peaks with fragmentation data',
	            checked: 'checked',
	            name: 'use_fragmentation'
	        },{
	            xtype: 'numberfield',
	            name: 'ms_intensity_cutoff',
	            fieldLabel: 'Absolute intensity minimum of lvl1 scan peaks which are matched with metabolites',
	            value: 200000.0,
	            decimalPrecision: 5
	        },{
	            xtype: 'numberfield',
	            name: 'msms_intensity_cutoff',
	            fieldLabel: 'Ratio of basepeak intensity',
	            value: 0.1,
	            decimalPrecision: 5
	        },{
	            xtype: 'numberfield',
	            name: 'mz_precision',
	            fieldLabel: 'M/z offset which is allowed for matching a metabolite mass to m/z of a peak',
	            value: 0.01,
	            decimalPrecision: 5
	        }]
	    }]
    }],
    buttons: [{
      text: 'Submit',
      handler: function(){
          var form = this.up('form').getForm();
          if(form.isValid()){
              form.submit({
                  url: '${request.route_url('home')}',
                  waitMsg: 'Uploading your data...',
                  success: function(fp, o) {
                      window.location = '${request.application_url}/status/'+o.result.jobid;
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
});
</script>
</head>
<body>
</body>
</html>