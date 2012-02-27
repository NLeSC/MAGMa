<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <title>MAGMa - Homepage</title>
<link rel="stylesheet" href="${request.extjsroot}/resources/css/ext-all.css" type="text/css"></link>
<link rel="stylesheet" href="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb.css')}" type="text/css"></link>
<link rel="stylesheet" href="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.8.7.custom.css')}" type="text/css"></link>
<script type="text/javascript" src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb-libs.js')}"></script>
<script type="text/javascript" src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb.js')}"></script>
<script type="text/javascript" src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.8.7.custom.min.js')}"></script>
<script type="text/javascript" src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/ChemDoodleWeb-sketcher.js')}"></script>
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
    width: 800,
    title: 'MAGMa submit form',
    bodyPadding: 5,
    renderTo: Ext.getBody(),
    defaults: { bodyPadding: 5 },
    items:[{
       xtype: 'fieldset',
       title: "Structures",
       items:[{
           xtype: 'tabpanel',
           activeTab: 0,
           defaults: { bodyPadding: 5 },
           plain: true,
           items: [{
               title: 'Upload',
               items:[{
                   xtype: 'combo',
                   id: 'structures_format',
                   store: ['smile', 'SD' ],
                   allowBlank: false,
                   value: 'smile',
               }, {
                   xtype: 'textarea',
                   name: 'structures',
                   id: 'structures_area',
                   emptyText: 'Enter smile string followed by space and name on each line',
                   height: 200,
                   width: 500,
               }, {
                   xtype: 'displayfield',
                   value: 'or'
               }, {
                   xtype: 'filefield',
                   emptyText: 'Upload structures in file',
                   width: 300
               }]
           }, {
               title: 'Draw',
               items: [{
                   id: 'sketcher', // during form submit fetch molblock from sketcher
                   contentEl: 'sketcher_content',
                   width: 500,
                   height: 300
               }]
           }]
       }],
    }, {
        xtype: 'fieldset',
        title: 'MS/MS data',
        layout: 'hbox',
        items: [{
            xtype: 'combo',
            store: ['mzxml', 'peaklist' ],
            allowBlank: false,
            value: 'mzxml'
        }, {
            name: 'db',
            xtype: 'filefield',
            emptyText: 'Upload MS/MS data file',
            width: 300
        }]
    }, {
        xtype: 'fieldset',
        title: 'Generate metabolite options',
        collapsed: true,
        collapsible: true,
        defaults: { labelWidth: 300 },
        items: [{
            fieldLabel: 'Maximum number of reaction steps',
            name: 'n_reaction_steps',
            xtype: 'numberfield',
            allowBlank: false,
            value: 2,
            maxValue: 10,
            minValue: 0,
            decimalPrecision: 0
        },{
            xtype:'combobox',
            fieldLabel: 'Metabolism types',
            store:['phase1', 'phase2'],
            multiSelect: true,
            allowBlank: false,
            value: 'phase1, phase2',
            name: 'metabolism_types'
        }]
      },{
        xtype: 'fieldset',
        title: 'MS data options',
        collapsed: true,
        collapsible: true,
        defaults: { labelWidth: 300 },
        items: [{
            xtype: 'numberfield',
            name: 'max_ms_level',
            fieldLabel: 'Maximum MS level',
            allowBlank: false,
            value: 3,
            minValue: 1,
            maxValue: 5,
            decimalPrecision: 0
        },{
            xtype: 'numberfield',
            name: 'precursor_mz_precision',
            fieldLabel: 'Precision for matching precursor mz with peak mz in parent scan',
            allowBlank: false,
            value: 0.001,
            decimalPrecision: 5
        },{
            xtype: 'numberfield',
            name: 'abs_peak_cutoff',
            fieldLabel: 'Absolute intensity threshold for storing peaks in database',
            allowBlank: false,
            value: 1000,
            decimalPrecision: 5
        },{
            xtype: 'numberfield',
            name: 'rel_peak_cutoff',
            fieldLabel: 'Fraction of basepeak intensity threshold threshold for storing peaks in database',
            allowBlank: false,
            value: 0.01,
            decimalPrecision: 5
        }]
    }, {
      xtype: 'fieldset',
      title: 'Annotate options',
      collapsed: true,
      collapsible: true,
      defaults: { labelWidth: 300 },
      items: [{
          fieldLabel: 'Maximum number of bonds broken in substructures generated from metabolites',
          name: 'max_broken_bonds',
          xtype: 'numberfield',
          allowBlank: false,
          value: 4,
          maxValue: 10,
          minValue: 0,
          decimalPrecision: 0
      },{
          fieldLabel: 'Ionisation mode',
          xtype: 'radiogroup',
          columns: 2,
          items: [
            { boxLabel:'Negative', checked: true, name:'ionisation_mode', inputValue: -1 },
            { boxLabel:'Positve', checked: false, name:'ionisation_mode', inputValue: 1 }
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
          name: 'use_msms_only'
      },{
          xtype: 'numberfield',
          name: 'ms_intensity_cutoff',
          fieldLabel: 'Absolute intensity minimum of lvl1 scan peaks which are matched with metabolites',
          allowBlank: false,
          value: 200000.0,
          decimalPrecision: 5
      },{
          xtype: 'numberfield',
          name: 'msms_intensity_cutoff',
          fieldLabel: 'Ratio of basepeak intensity',
          allowBlank: false,
          value: 0.1,
          decimalPrecision: 5
      },{
          xtype: 'numberfield',
          name: 'mz_precision',
          fieldLabel: 'M/z offset which is allowed for matching a metabolite mass to m/z of a peak',
          allowBlank: false,
          value: 0.01,
          decimalPrecision: 5
      }]
    }],
    buttons: [{
      text: 'Submit',
      handler: function(){
          var form = this.up('form').getForm();
          var mol = sketcher.getMolecule();
          if (mol.bonds.length > 0) {
              form.setValues({
                 structures_format: 'SD',
                 structures_area: ChemDoodle.writeMOL(mol)
              });
          }
          if(form.isValid()){
              form.submit({
                  url: '${request.route_url('home')}',
                  waitMsg: 'Uploading your data...',
                  success: function(fp, o) {
                      window.location = '${request.application_url}/status/'+o.result.jobid;
                  },
                  failure: function(form, action) {
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
});
</script>
</head>
<body>
<div id="sketcher_content" class="x-hidden">
<script language="javascript">
var sketcher = new ChemDoodle.SketcherCanvas(
        'sketcher_canvas', 500, 300,
        '${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/icons/')}',
        ChemDoodle.featureDetection.supports_touch(), false);
sketcher.repaint();
sketcher.toolbarManager.buttonSave.disable();
sketcher.toolbarManager.buttonOpen.disable();
</script>
</div>
</body>
</html>
