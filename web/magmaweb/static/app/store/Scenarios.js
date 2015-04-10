/**
 * Store for metabolize scenarios.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.store.Scenarios', {
  extend: 'Ext.data.Store',
  storeId: 'scenarios',
  fields: ['type', 'steps'],
  proxy: {
    type: 'memory',
    reader: 'json'
  }
});