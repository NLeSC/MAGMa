/**
 * Panel for metabolites contains grid and add form.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.metabolite.Panel', {
  extend: 'Ext.panel.Panel',
  alias: 'widget.metabolitepanel',
  id: 'metabolitepanel',
  title: 'Molecules',
  requires: [
    'Esc.magmaweb.view.metabolite.List',
    'Esc.magmaweb.view.metabolite.AddForm'
  ],
  tools: [{
     type: 'save',
     tooltip: 'Save metabolites as comma seperated file',
     action: 'download'
  }, {
     type: 'gear',
     tooltip: 'Perform actions on metabolites',
     action: 'actions'
  }],
  layout: 'card',
  items: [{
      xtype: 'metabolitelist'
  }, {
      xtype: 'metaboliteaddform',
      bodyPadding: 10
  }],
  /**
   * Shortcut to {@link Ext.layout.container.Card#setActiveItem setActiveItem} in layout.
   *
   * @param {Number} newCard 0 == grid, 1 == add form
   */
  setActiveItem: function(newCard) {
      this.getLayout().setActiveItem(newCard);
  }
});
