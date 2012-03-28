Ext.define('Esc.magmaweb.view.fragment.AnnotateFieldSet', {
    extend : 'Ext.form.FieldSet',
    alias : 'widget.annotatefieldset',
    requires : ['Ext.form.field.Number', 'Ext.form.RadioGroup',
            'Ext.form.field.Checkbox'],
    title : 'Annotate options',
    defaults : {
        labelWidth : 300
    },
    items : [{
        fieldLabel : 'Maximum number of bonds broken in substructures generated from metabolites',
        name : 'max_broken_bonds',
        xtype : 'numberfield',
        allowBlank : false,
        value : 4,
        maxValue : 10,
        minValue : 0,
        decimalPrecision : 0
    }, {
        fieldLabel : 'Ionisation mode',
        xtype : 'radiogroup',
        columns : 2,
        items : [{
                    boxLabel : 'Negative',
                    checked : true,
                    name : 'ionisation_mode',
                    inputValue : -1
                }, {
                    boxLabel : 'Positve',
                    checked : false,
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
    }, {
        xtype : 'numberfield',
        name : 'ms_intensity_cutoff',
        fieldLabel : 'Absolute intensity minimum of lvl1 scan peaks which are matched with metabolites',
        allowBlank : false,
        value : 200000.0,
        decimalPrecision : 5
    }, {
        xtype : 'numberfield',
        name : 'msms_intensity_cutoff',
        fieldLabel : 'Ratio of basepeak intensity',
        allowBlank : false,
        value : 0.1,
        decimalPrecision : 5
    }, {
        xtype : 'numberfield',
        name : 'mz_precision',
        fieldLabel : 'M/z offset which is allowed for matching a metabolite mass to m/z of a peak',
        allowBlank : false,
        value : 0.001,
        decimalPrecision : 5
    }]
});