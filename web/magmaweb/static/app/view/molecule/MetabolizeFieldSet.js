/**
 * Fieldset with metabolize options
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.molecule.MetabolizeFieldSet', {
  extend: 'Ext.form.FieldSet',
  alias: 'widget.metabolizefieldset',
  requires: [
    'Ext.form.field.Checkbox',
    'Ext.form.FieldContainer',
    'Esc.magmaweb.view.molecule.ScenarioField'
  ],
  title: 'Metabolize',
  checkboxName: 'metabolize',
  checkboxToggle: true,
  frame: true,
  items: [{
    xtype: 'displayfield',
    hidden: true,
    itemId: 'onemol_metabolize_warning',
    value: 'Only first molecule will be metabolized'
  }, {
    xtype: 'fieldcontainer',
    items: [{
      title: 'Scenario',
      xtype: 'scenariofield',
      name: 'scenario'
    }]
  }]
});
