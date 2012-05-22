/**
 * Metabolite model.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.model.Metabolite', {
  idProperty: 'metid',
  extend:'Ext.data.Model',
  fields: [{
    name: 'metid'
  },{
    name: 'mol', type:'string'
  },{
    name: 'level'
  },{
    name: 'probability'
  },{
    name: 'reactionsequence'
  },{
    name: 'smiles'
  },{
    name: 'molformula'
  },{
    name: 'isquery', type: 'bool'
  },{
    name: 'origin'
  },{
    name: 'nhits'
  },{
    name: 'atoms', defaultValue: [] // array of atom indexes of molecule which are the substructure of the query
  },{
    name: 'nr_scans', type:'number'
  },{
    name: 'scans', defaultValue: [] // Filled when metabolite is selected
  },{
    name: 'score', type:'number' // Only filled when scan is selected
  }, {
    name: 'mim', type:'number'
  }, {
    name: 'logp', type:'number'
  }]
});
