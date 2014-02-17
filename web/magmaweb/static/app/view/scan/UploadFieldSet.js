/**
 * Fieldset with ms data upload options.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.scan.UploadFieldSet', {
    extend : 'Ext.form.Panel',
    alias: 'widget.uploadmsdatafieldset',
    requires: [
         'Esc.form.field.TextareaTab',
         'Ext.form.field.ComboBox',
         'Ext.form.field.Number',
         'Ext.form.field.File',
         'Ext.container.Container'
    ],
    frame: true,
    bodyPadding: '5',
    layout: 'anchor',
    defaults : {
        labelWidth : 120
    },
    items: [{
        xtype: 'combo',
        store: [['mzxml','mzXML'],
                ['mass_tree', 'Mass Tree'],
                ['form_tree', 'Formula Tree']],
        allowBlank: false,
        fieldLabel: 'Format',
        name: 'ms_data_format',
        value: 'mzxml',
        anchor: '75%'
    }, {
        xtype : 'textareatab',
        name : 'ms_data',
        id: 'ms_data_area',
        emptyText : 'Enter MS data in a Tree format or mzXML',
        height : 200,
        anchor : '90%'
    }, {
        xtype: 'container',
        layout: 'hbox',
        items: [{
            xtype: 'displayfield',
            value : '<a title="Information on examples">Examples:</a>'
        }, {
            xtype: 'button',
            margin: '0 5 0 5',
            text: 'Chlorogenic acid (Mass Tree)',
            flex: 1,
            action: 'loadmsdataexample'
        }, {
            xtype: 'button',
            text: 'Chlorogenic acid (Formula Tree)',
            flex: 1,
            action: 'loadmsdataexample2'
        }],
        anchor : '90%'
    }, {
        xtype : 'displayfield',
        value : '<br>or'
    }, {
        name: 'ms_data_file',
        xtype: 'filefield',
        emptyText: 'Upload MS/MS data file',
        anchor: '75%'
    }, {
        xtype: 'displayfield',
        name: 'filter_heading',
        value: '<br>Filter options:'
    }, {
        xtype: 'numberfield',
        name: 'max_ms_level',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&le;</span>',
        fieldLabel: 'MS level',
        labelAttrTpl: 'data-qwidth=100 data-qtip="Maximum MS level"',
        allowBlank: false,
        minValue: 1,
        maxValue: 15,
        allowDecimals: false,
        anchor: '75%'
    }, {
        xtype: 'numberfield',
        name: 'abs_peak_cutoff',
        fieldLabel: 'Noise filter',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&lt;</span>',
        labelAttrTpl: 'data-qwidth=200 data-qtip="Absolute intensity threshold for storing peaks in database"',
        allowBlank: false,
        decimalPrecision: 5,
        anchor: '75%'
    }, {
        xtype: 'numberfield',
        fieldLabel: 'MS1 scan number',
        labelSeparator: '',
        labelAttrTpl: 'data-qwidth=200 data-qtip="Read only spectral tree specified by MS1 scan number"',
        name: 'scan',
        allowBlank: true,
        minValue: 0,
        allowDecimals: false,
        anchor: '75%'
    }]
});
