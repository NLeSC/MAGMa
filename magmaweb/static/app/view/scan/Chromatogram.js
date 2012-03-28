/**
 * Chromatogram of lc-ms.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.scan.Chromatogram', {
  extend: 'Ext.panel.Panel',
  alias: 'widget.chromatogrampanel',
  title: 'Chromatogram',
  requires: ['Esc.d3.Chromatogram'],
  tools:[{
    type:'search',
    tooltip: 'Select lvl1 scan by identifier',
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
    tooltip: 'Upload MS data',
    action: 'upload'
  }],
  layout: 'fit',
  items: {
      emptyText: 'No chromatogram available: Upload ms data',
      axesPadding: [16, 5, 58, 80],
      xtype: 'chromatogram'
  }
});
