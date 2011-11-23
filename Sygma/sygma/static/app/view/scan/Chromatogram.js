Ext.define('Esc.msygma.view.scan.Chromatogram', {
  extend: 'Ext.esc.Chromatogram',
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
