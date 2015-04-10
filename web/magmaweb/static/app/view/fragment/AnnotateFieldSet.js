/**
 * Fieldset with annotate options.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.fragment.AnnotateFieldSet', {
  extend: 'Ext.form.Panel',
  alias: 'widget.annotatefieldset',
  requires: [
    'Ext.form.field.Number', 'Ext.form.RadioGroup',
    'Ext.form.field.Checkbox', 'Ext.form.field.Radio'
  ],
  title: 'Annotate',
  frame: true,
  bodyPadding: '5',
  layout: 'anchor',
  defaults: {
    labelWidth: 150
  },
  items: [{
    fieldLabel: 'Ionisation mode',
    xtype: 'radiogroup',
    columns: 1,
    vertical: true,
    items: [{
      boxLabel: 'Negative',
      name: 'ionisation_mode',
      inputValue: -1
    }, {
      boxLabel: 'Positve',
      name: 'ionisation_mode',
      inputValue: 1
    }],
    anchor: '75%'
  }, {
    xtype: 'displayfield',
    value: '<br>Substructure options:'
  }, {
    fieldLabel: 'Bond dissociations',
    labelAttrTpl: 'data-qwidth=200 data-qtip="Maximum number of bond breaks to generate substructures"',
    labelSeparator: '',
    afterLabelTextTpl: '<span class="relation">&le;</span>',
    name: 'max_broken_bonds',
    xtype: 'numberfield',
    allowBlank: false,
    maxValue: 4,
    minValue: 0,
    decimalPrecision: 0,
    anchor: '75%'
  }, {
    fieldLabel: 'Additional small losses',
    labelAttrTpl: 'data-qwidth=200 data-qtip="Maximum number of additional water (OH) and/or ammonia (NH2) losses"',
    labelSeparator: '',
    afterLabelTextTpl: '<span class="relation">&le;</span>',
    name: 'max_water_losses',
    xtype: 'numberfield',
    allowBlank: false,
    maxValue: 4,
    minValue: 0,
    decimalPrecision: 0,
    anchor: '75%'
  }, {
    xtype: 'displayfield',
    name: 'precision_heading',
    value: '<br>Accuracy:'
  }, {
    xtype: 'numberfield',
    name: 'mz_precision',
    labelSeparator: '',
    afterLabelTextTpl: '<span class="relation">&le;</span>',
    fieldLabel: 'Relative (ppm)',
    labelAttrTpl: 'data-qwidth=200 data-qtip="Maximum relative m/z error (ppm)"',
    allowBlank: false,
    maxValue: 1000,
    minValue: 0,
    decimalPrecision: 5,
    anchor: '75%'
  }, {
    xtype: 'numberfield',
    name: 'mz_precision_abs',
    labelSeparator: '',
    afterLabelTextTpl: '<span class="relation">&le;</span>',
    fieldLabel: 'Absolute (Da)',
    labelAttrTpl: 'data-qwidth=200 data-qtip="Maximum absolute m/z error (Da)"',
    allowBlank: false,
    maxValue: 1,
    minValue: 0,
    decimalPrecision: 5,
    anchor: '75%'
  }, {
    xtype: 'numberfield',
    name: 'precursor_mz_precision',
    labelSeparator: '',
    afterLabelTextTpl: '<span class="relation">&le;</span>',
    fieldLabel: 'Precursor m/z (Da)',
    labelAttrTpl: 'data-qwidth=200 data-qtip="Maximum absolute error of precursor m/z values"',
    allowBlank: false,
    maxValue: 1,
    minValue: 0,
    decimalPrecision: 5,
    anchor: '75%'
  }, {
    xtype: 'displayfield',
    name: 'intensity_heading',
    value: '<br>Intensity thresholds:'
  }, {
    xtype: 'numberfield',
    name: 'ms_intensity_cutoff',
    labelSeparator: '',
    afterLabelTextTpl: '<span class="relation">&ge;</span>',
    fieldLabel: 'MS<sup>1</sup> (abs.)',
    labelAttrTpl: 'data-qwidth=200 data-qtip="Minimum intensity of MS1 precursor ion peaks to be annotated"',
    allowBlank: false,
    decimalPrecision: 5,
    anchor: '75%'
  }, {
    xtype: 'numberfield',
    name: 'msms_intensity_cutoff',
    labelSeparator: '',
    afterLabelTextTpl: '<span class="relation">&ge;</span>',
    fieldLabel: 'MS<sup>n&gt;1</sup> (% of base peak)',
    labelAttrTpl: 'data-qwidth=200 data-qtip="Minimum intensity of fragment peaks to be annotated, as percentage of basepeak"',
    allowBlank: false,
    decimalPrecision: 5,
    anchor: '75%'
  }]
});