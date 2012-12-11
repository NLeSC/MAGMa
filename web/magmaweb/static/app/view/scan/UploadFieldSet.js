/**
 * Fieldset with ms data upload options.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.scan.UploadFieldSet', {
	extend : 'Ext.form.FieldSet',
    alias: 'widget.uploadmsdatafieldset',
    requires: [
         'Ext.form.field.ComboBox',
         'Ext.form.field.Number',
         'Ext.form.field.File',
         'Ext.container.Container'
    ],
    items: [{
        xtype: 'combo',
        store: [['mzxml','mzXML'], ['tree', 'Tree']],
        allowBlank: false,
        fieldLabel: 'Format',
        name: 'ms_data_format',
        value: 'mzxml'
    }, {
        xtype : 'textarea',
        name : 'ms_data',
        id: 'ms_data_area',
        emptyText : 'Enter MS data in Tree format or mzXML',
        height : 200,
        width : 500
    }, {
        xtype : 'displayfield',
        value : 'or'
    }, {
        name: 'ms_data_file',
        xtype: 'filefield',
        emptyText: 'Upload MS/MS data file',
        width: 300
    }, {
        xtype : 'displayfield',
        value : 'or'
    }, {
    	xtype: 'container',
    	layout: 'hbox',
    	items: [{
        	xtype: 'button',
	    	text: 'Load chlorogenic acid example',
			tooltip: 'Loads chlorogenic acid example ms data set and configuration which gives well annotated result',
	    	action: 'loadmsdataexample'
	    }, {
            xtype : 'displayfield',
	        flex: 1,
	        value : '&nbsp;<a href="help/example">Information on example</a>'
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
            minValue: 1,
            maxValue: 15,
            decimalPrecision: 0
        }, {
            xtype: 'numberfield',
            name: 'abs_peak_cutoff',
            fieldLabel: 'Absolute intensity threshold for storing peaks in database',
            allowBlank: false,
            decimalPrecision: 5
        }, {
            xtype: 'numberfield',
            name: 'rel_peak_cutoff',
            fieldLabel: 'Fraction of basepeak intensity threshold for storing peaks in database',
            allowBlank: false,
            decimalPrecision: 5
        }]
    }]
});
