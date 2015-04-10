/**
 * Form to add molecules.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.molecule.AddForm', {
  extend: 'Ext.form.Panel',
  alias: 'widget.moleculeaddform',
  autoScroll: true,
  trackResetOnLoad: true,
  requires: [
    'Esc.magmaweb.view.molecule.AddFieldSet',
    'Esc.magmaweb.view.molecule.MetabolizeFieldSet',
    'Esc.magmaweb.view.fragment.AnnotateFieldSet'
  ],
  items: [{
    margin: '0 0 10 0',
    xtype: 'addstructurefieldset'
  }, {
    xtype: 'metabolizefieldset',
    margin: '0 0 10 0',
    checkboxToggle: true,
    checkboxName: 'metabolize',
    collapsed: true,
    collapsible: true
  }, {
    xtype: 'annotatefieldset',
    collapsed: true,
    collapsible: true
  }],
  buttons: [{
    text: 'Add structures',
    action: 'addstructures'
  }, {
    text: 'Reset',
    handler: function() {
      this.up('form').getForm().reset();
    }
  }, {
    text: 'Cancel',
    action: 'addstructurescancel'
  }],
  /**
   * Enable or disable annotate fieldset.
   *
   * @param {Boolean} disabled True to disable
   */
  setDisabledAnnotateFieldset: function(disabled) {
    this.query('annotatefieldset')[0].setDisabled(disabled);
    this.queryById('structure_database_tab').setDisabled(disabled);
  },
  /**
   * Load form defaults from server.
   *
   * @param {String} url Url to load defaults from.
   */
  loadDefaults: function(url) {
    this.load({
      url: url,
      method: 'GET',
      waitMsg: 'Fetching defaults',
      failure: function(form, action) {
        Ext.Error.raise(action.response.responseText);
      }
    });
  }
});
