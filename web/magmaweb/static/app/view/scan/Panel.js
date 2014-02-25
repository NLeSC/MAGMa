/**
 * Chromatogram of lc-ms and ms data upload form in card layout.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.scan.Panel', {
  extend: 'Ext.panel.Panel',
  alias: 'widget.scanpanel',
  title: 'Chromatogram',
  id: 'scanpanel',
  requires: [
    'Esc.d3.Chromatogram',
    'Esc.magmaweb.view.scan.UploadForm'
  ],
  tools:[{
    type:'search',
    tooltip: 'Select level 1 scan by identifier',
    action: 'search'
  }, {
    type:'refresh',
    tooltip:'Clear scan selection',
    action: 'clearselection'
  }, {
    type: 'restore',
    tooltip: 'Center chromatogram',
    action: 'center'
  }, {
    type: 'gear',
    tooltip: 'Upload MS data & Change zoom direction',
    action: 'actions'
  }, {
    type: 'help',
    tooltip: 'Help',
    action: 'help'
  }],
  layout: 'card',
  items: [{
      emptyText: 'No chromatogram available: Upload ms data',
      axesPadding: [16, 5, 58, 80],
      xtype: 'chromatogram'
  }, {
      xtype: 'scanuploadform',
      bodyPadding: 10
  }],
  /**
   * Shortcut to {@link Ext.layout.container.Card#setActiveItem setActiveItem} in layout.
   *
   * @param {Number} newCard 0 == chromatogram, 1 == upload form
   */
  setActiveItem: function(newCard) {
      this.getLayout().setActiveItem(newCard);
  }
});
