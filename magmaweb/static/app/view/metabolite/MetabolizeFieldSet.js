/**
 * Fieldset with metabolize options
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.metabolite.MetabolizeFieldSet', {
    extend : 'Ext.form.FieldSet',
    alias : 'widget.metabolizefieldset',
    requires : [ 'Ext.form.field.Number', 'Ext.form.field.ComboBox' ],
    title : 'Generate metabolite options',
    defaults : {
        labelWidth : 300
    },
    items : [ {
        fieldLabel : 'Maximum number of reaction steps',
        name : 'n_reaction_steps',
        xtype : 'numberfield',
        allowBlank : false,
        value : 2,
        maxValue : 10,
        minValue : 0,
        decimalPrecision : 0
    }, {
        xtype : 'combobox',
        fieldLabel : 'Metabolism types',
        store : [ 'phase1', 'phase2' ],
        multiSelect : true,
        allowBlank : false,
        value : [ 'phase1', 'phase2' ],
        name : 'metabolism_types'
    } ]
});