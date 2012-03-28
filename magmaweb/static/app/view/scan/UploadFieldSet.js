Ext.define('Esc.magmaweb.view.scan.UploadFieldSet', {
    extend: 'Ext.container.Container',
    alias: 'widget.uploadmsdatafieldset',
    requires: [
         'Ext.form.field.ComboBox',
         'Ext.form.field.Number',
         'Ext.form.field.File',
         'Ext.form.FieldSet'
    ],
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
    }]
});