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
        fieldLabel : 'Maximum number of bonds broken in substructures generated from metabolites',
        name : 'max_broken_bonds',
        xtype : 'numberfield',
        allowBlank : false,
        maxValue : 10,
        minValue : 0,
        decimalPrecision : 0
    }, {
        fieldLabel : 'Ionisation mode',
        xtype : 'radiogroup',
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
        xtype : 'checkbox',
        fieldLabel : 'Skip fragmentation',
        name : 'skip_fragmentation'
    }, {
        xtype : 'checkbox',
        fieldLabel : 'Annotate all lvl1 peaks, including those without fragmentation data',
        name : 'use_all_peaks'
    },{
        xtype: 'numberfield',
        name: 'precursor_mz_precision',
        fieldLabel: 'Precision for matching precursor mz with peak mz in parent scan',
        allowBlank: false,
        decimalPrecision: 5
    }, {
        xtype : 'numberfield',
        name : 'ms_intensity_cutoff',
        fieldLabel : 'Absolute intensity minimum of lvl1 scan peaks which are matched with metabolites',
        allowBlank : false,
        decimalPrecision : 5
    }, {
        xtype : 'numberfield',
        name : 'msms_intensity_cutoff',
        fieldLabel : 'Ratio of basepeak intensity',
        allowBlank : false,
        decimalPrecision : 5
    }, {
        xtype : 'numberfield',
        name : 'mz_precision',
        fieldLabel : 'M/z offset which is allowed for matching a metabolite mass to m/z of a peak',
        allowBlank : false,
        decimalPrecision : 5
    }]
});