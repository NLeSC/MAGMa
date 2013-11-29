/**
 * Fieldset with metabolize options
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.metabolite.MetabolizeFieldSet', {
    extend : 'Ext.panel.Panel',
    alias : 'widget.metabolizefieldset',
    requires : [ 'Ext.form.field.Number', 'Ext.form.field.ComboBox' ],
    title : 'Metabolize',
    frame: true,
    defaults : {
        labelWidth : 300
    },
    items : [ {
        fieldLabel : 'Perform metabolization of molecules',
        xtype: 'checkbox',
        name: 'metabolize'
    }, {
        fieldLabel : 'Maximum number of reaction steps',
        name : 'n_reaction_steps',
        xtype : 'numberfield',
        allowBlank : false,
        maxValue : 10,
        minValue : 0,
        decimalPrecision : 0
    }, {
        xtype : 'combobox',
        fieldLabel : 'Metabolism types',
        store : [ 'digest', 'gut', 'phase1', 'phase2' ],
        multiSelect : true,
        allowBlank : false,
        name : 'metabolism_types'
    } ]
});