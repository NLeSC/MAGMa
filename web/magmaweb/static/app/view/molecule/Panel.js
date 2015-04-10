/**
 * Panel for molecules contains grid and add form.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.molecule.Panel', {
  extend: 'Ext.panel.Panel',
  alias: 'widget.moleculepanel',
  id: 'moleculepanel',
  title: 'Molecules',
  requires: [
    'Esc.magmaweb.view.molecule.List',
    'Esc.magmaweb.view.molecule.AddForm'
  ],
  tools: [{
    type: 'save',
    tooltip: 'Save molecules',
    action: 'download'
  }, {
    type: 'gear',
    tooltip: 'Perform actions on molecules',
    action: 'actions'
  }, {
    type: 'help',
    tooltip: 'Help',
    action: 'help'
  }],
  layout: 'card',
  items: [{
    xtype: 'moleculelist'
  }, {
    xtype: 'moleculeaddform',
    bodyPadding: 10
  }],
  /**
   * Shortcut to {@link Ext.layout.container.Card#setActiveItem setActiveItem} in layout.
   *
   * @param {Number} newCard 0 == grid, 1 == add form
   */
  setActiveItem: function(newCard) {
    this.getLayout().setActiveItem(newCard);
  },
  getPagingToolbar: function() {
    return this.items.getAt(0).dockedItems.getAt(1);
  }
});
