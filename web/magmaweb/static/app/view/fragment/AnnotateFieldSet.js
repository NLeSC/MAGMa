/**
 * Fieldset with annotate options.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.fragment.AnnotateFieldSet', {
    extend : 'Ext.form.FieldSet',
    alias : 'widget.annotatefieldset',
    requires : [
        'Ext.form.field.Number', 'Ext.form.RadioGroup',
        'Ext.form.field.Checkbox', 'Ext.form.field.Radio'
    ],
    title : 'Annotate',
    defaults : {
        labelWidth : 200
    },
    items : [{
        fieldLabel : 'Ionisation mode',
        xtype : 'radiogroup',
        width: 500,
        columns : 2,
        items : [{
                    boxLabel : 'Negative',
                    name : 'ionisation_mode',
                    inputValue : -1
                }, {
                    boxLabel : 'Positve',
                    name : 'ionisation_mode',
                    inputValue : 1
                }]
    }, {
        xtype: 'displayfield',
        value: 'Substructure options:'
    }, {
        fieldLabel: 'Bond breaks',
        tooltip : 'Maximum number of bond breaks to generate substructures',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&le;</span>',
        name : 'max_broken_bonds',
        xtype : 'numberfield',
        allowBlank : false,
        maxValue : 10,
        minValue : 0,
        decimalPrecision : 0
    }, {
        fieldLabel: 'Additional water losses',
        tooltip : 'Maximum number of additional neutral water losses',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&le;</span>',
        name : 'max_water_losses',
        xtype : 'numberfield',
        allowBlank : false,
        maxValue : 10,
        minValue : 0,
        decimalPrecision : 0
    }, {
        xtype: 'displayfield',
        value: 'Precision:'
    }, {
        xtype : 'numberfield',
        name : 'mz_precision',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&le;</span>',
        fieldLabel: 'Relative (ppm)',
        tooltip : 'Mass precision (ppm) for matching calculated masses with peaks',
        allowBlank : false,
        decimalPrecision : 5
    }, {
        xtype : 'numberfield',
        name : 'mz_precision_abs',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&le;</span>',
        fieldLabel: 'Absolute (Da)',
        tooltip : 'Mass precision (Da) for matching calculated masses with peaks',
        allowBlank : false,
        decimalPrecision : 5
    },{
        xtype: 'numberfield',
        name: 'precursor_mz_precision',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&le;</span>',
        fieldLabel: 'Precursor m/z (Da)',
        tooltip: 'Mass precision for matching peaks and precursor ions',
        allowBlank: false,
        decimalPrecision: 5
    }, {
        xtype: 'displayfield',
        value: 'Intensity thresholds:'
    }, {
        xtype : 'numberfield',
        name : 'ms_intensity_cutoff',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&le;</span>',
        fieldLabel: 'MS<sup>1</sup> (abs.)',
        tooltip : 'Minimum intensity of level 1 peaks to be annotated',
        allowBlank : false,
        decimalPrecision : 5
    }, {
        xtype : 'numberfield',
        name : 'msms_intensity_cutoff',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&le;</span>',
        fieldLabel: 'MS<sup>n&gt;1</sup> (% of base peak)',
        tooltip : 'Minimum intensity of fragment peaks to be annotated, as percentage of basepeak',
        allowBlank : false,
        decimalPrecision : 5
    }]
});
