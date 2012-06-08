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
    title : 'Annotate options',
    defaults : {
        labelWidth : 300
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
    },{
        fieldLabel : 'Maximum number of bond breaks to generate substructures',
        name : 'max_broken_bonds',
        xtype : 'numberfield',
        allowBlank : false,
        maxValue : 10,
        minValue : 0,
        decimalPrecision : 0
    }, {
        xtype : 'numberfield',
        name : 'mz_precision',
        fieldLabel : 'Mass precision for matching calculated masses with peaks',
        allowBlank : false,
        decimalPrecision : 5
    },{
        xtype: 'numberfield',
        name: 'precursor_mz_precision',
        fieldLabel: 'Mass precision for matching peaks and precursor ions',
        allowBlank: false,
        decimalPrecision: 5
    }, {
        xtype : 'numberfield',
        name : 'ms_intensity_cutoff',
        fieldLabel : 'Minimum intensity of level 1 peaks to be annotated',
        allowBlank : false,
        decimalPrecision : 5
    }, {
        xtype : 'numberfield',
        name : 'msms_intensity_cutoff',
        fieldLabel : 'Minimum intensity of fragment peaks to be annotated, as fraction of basepeak',
        allowBlank : false,
        decimalPrecision : 5
    }, {
        xtype : 'checkbox',
        fieldLabel : 'Annotate all level 1 peaks, including those not fragmented',
        name : 'use_all_peaks'
    },  {
        xtype : 'checkbox',
        fieldLabel : 'Skip fragmentation',
        name : 'skip_fragmentation'
    }]
});
