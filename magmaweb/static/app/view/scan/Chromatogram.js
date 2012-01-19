/**
 * Chromatogram of lc-ms.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.scan.Chromatogram', {
  extend: 'Esc.d3.Chromatogram',
  alias: 'widget.scanchromatogram',
  title: 'Chromatogram',
  emptyText: 'Loading chromatogram ...',
  axesPadding: [16, 5, 58, 80],
  tools:[{
    type:'search',
    tooltip: 'Select lvl1 scan by identifier',
    action: 'search'
  },{
    type:'refresh',
    tooltip:'Clear scan selection',
    action: 'clearselection',
  }]
});
