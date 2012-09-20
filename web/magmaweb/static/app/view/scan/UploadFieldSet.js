/**
 * Fieldset with ms data upload options.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
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
        store: [['mzxml','mzXML']],
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
            minValue: 1,
            maxValue: 15,
            decimalPrecision: 0
        },{
            xtype: 'numberfield',
            name: 'abs_peak_cutoff',
            fieldLabel: 'Absolute intensity threshold for storing peaks in database',
            allowBlank: false,
            decimalPrecision: 5
        }]
    }]
});
