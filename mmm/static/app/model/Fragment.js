/**
 * Fragment model.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.mmm.model.Fragment', {
  extend:'Ext.data.Model',
  idProperty: 'fragid',
  fields: [{
    name: 'fragid',
  },{
    name: 'metid',
  },{
    name: 'scanid'
  },{
    name: 'mz'
  },{
    name: 'mass'
  },{
    name: 'score'
  },{
    name: 'parentfragid'
  },{
    name: 'mol', type: 'string'
  },{
    name: 'atoms', defaultValue: [] // array of atom indexes of molecule which are the substructure of the query
  },{
    name: 'deltah'
  },{
    name: 'mslevel'
  }],
  hasMany: { model: 'Fragment', name:'children' }
});

