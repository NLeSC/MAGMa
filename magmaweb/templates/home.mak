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
<style type="text/css">
.x-logo a {
  font-size: 40px;
  font-weight: bold;
  color: #333;
  text-decoration:none;
  padding-left: 520px;
  padding-top: 3px; /* aligns app title with text in logo  */
  background: url(${request.static_url('magmaweb:static/ESCIENCE_log_B_nl_long_cyanblack.jpg')}) no-repeat 5px 4px;
}

#welcome h1 {
  font-size: 200%;
}

</style>
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
    border: false,
    region: 'center',
    bodyPadding: 5,
    defaults: { bodyPadding: 5 },
    autoScroll: true,
    items:[{
        contentEl: 'welcome',
        border: false
    },{
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
                   name: 'structure_format',
                   store: ['smiles', 'sdf' ],
                   allowBlank: false,
                   value: 'smiles',
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
                   name: 'structures_file',
                   id: 'structures_filefield',
                   emptyText: 'Upload structures in file',
                   width: 300
               }]
           }, {
               title: 'Draw',
               items: [{
                   id: 'sketcher', // during form submit fetch molblock from sketcher
                   contentEl: 'sketcher_content',
                   border: false
               }]
           }]
       }],
    }, {
        xtype: 'fieldset',
        title: 'MS/MS data',
        layout: 'hbox',
        items: [{
            xtype: 'combo',
            store: ['mzxml'],
            allowBlank: false,
            name: 'ms_data_format',
            value: 'mzxml'
        }, {
            name: 'ms_data_file',
            xtype: 'filefield',
            allowBlank: false,
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
          fieldLabel: 'Skip fragmentation',
          name: 'skip_fragmentation'
      },{
          xtype:'checkbox',
          fieldLabel: 'Annotate all lvl1 peaks, including those without fragmentation data',
          name: 'use_all_peaks'
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
          value: 0.001,
          decimalPrecision: 5
      }],
    }, {
        xtype: 'textfield',
        name: 'description',
        fieldLabel: 'Description',
        width: 510,
        emptyText: 'Enter optional description which can be used to identify the results'
    }],
    buttons: [{
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

  Ext.create('Ext.container.Viewport', {
    layout: 'border',
    items: [{
      border: false,
      region: 'north',
      layout: {
        type: 'hbox',
        align: 'middle',
        padding: 2
      },
      items: [{
        xtype: 'component',
        cls: 'x-logo',
        html: '<a href="${request.route_url('home')}" data-qtip="<b>M</b>s <b>A</b>nnotation based on in silico <b>G</b>enerated <b>M</b>et<b>a</b>bolites">MAGMa</a>'
      }, {
        xtype:'tbspacer',
        flex:1 // aligns buttongroup right
      }, {
          xtype: 'buttongroup',
          columns: 3,
          items: [{
            text: 'Help',
            tooltip: 'Goto help pages',
            disabled: true
          },{
            text: 'Upload result',
            tooltip: 'Upload a result db for viewing',
            url: '${request.route_url('uploaddb')}'
          }]
        }]
     },
     form
     ]
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
<div id="welcome" class="x-hidden">
<h1>Welcome to the <b>M</b>s <b>A</b>nnotation based on in silico <b>G</b>enerated <b>M</b>et<b>a</b>bolites application</h1>
</div>
</body>
</html>
