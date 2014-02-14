/**
 * Fieldset with metabolize options
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.metabolite.MetabolizeFieldSet', {
    extend: 'Ext.form.FieldSet',
    alias: 'widget.metabolizefieldset',
    requires: [
         'Ext.form.field.Checkbox',
         'Ext.form.FieldContainer',
         'Esc.magmaweb.view.metabolite.ScenarioField'
    ],
    title: 'Metabolize',
    frame: true,
    defaults: {
        labelWidth: 200
    },
    items: [ {
        fieldLabel: 'Perform metabolization of molecules',
        xtype: 'checkbox',
        name: 'metabolize',
        boxLabel: '&nbsp;'
    }, {
        xtype: 'fieldcontainer',
        items: [{
            title: 'Scenario',
            xtype: 'scenariofield',
            name: 'scenario'
        }]
    }]
});