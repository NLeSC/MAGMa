/**
 * Windowed form to annotate structures and ms data.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.fragment.AnnotateForm', {
  extend: 'Ext.window.Window',
  title: 'Annotate all structures and ms data',
  alias: 'widget.annotateform',
  requires: [
    'Ext.form.Panel',
    'Esc.magmaweb.view.fragment.AnnotateFieldSet'
  ],
  modal: true,
  height: 400,
  width: 600,
  layout: 'fit',
  closeAction: 'hide',
  items: {
    xtype: 'form',
    trackResetOnLoad: true,
    bodyPadding: 5,
    defaults: {
      bodyPadding: 5
    },
    border: false,
    autoScroll: true,
    items: [{
      xtype: 'annotatefieldset'
    }],
    buttons: [{
      text: 'Annotate',
      action: 'annotate'
    }, {
      text: 'Reset',
      handler: function() {
        this.up('form').getForm().reset();
      }
    }]
  },
  /**
   * Get form.
   *
   * @return {Ext.form.Basic}
   */
  getForm: function() {
    return this.down('form').getForm();
  },
  /**
   * Load form defaults from server.
   * During loading this.loading will be true.
   *
   * @param {String} url Url to load defaults from.
   */
  loadDefaults: function(url) {
    var me = this;
    this.loading = true;
    this.getForm().load({
      url: url,
      method: 'GET',
      waitMsg: 'Fetching defaults',
      failure: function(form, action) {
        delete me.loading;
        Ext.Error.raise(action.response.responseText);
      },
      success: function(form) {
        delete me.loading;
      }
    });
  }
});