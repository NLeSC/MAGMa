/**
 * Molecule model.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.model.Molecule', {
  idProperty: 'molid',
  extend:'Ext.data.Model',
  fields: [{
    name: 'molid'
  },{
    name: 'mol', type:'string'
  },{
    name: 'level'
  },{
    name: 'refscore'
  },{
    name: 'reactionsequence'
  },{
    name: 'inchikey14'
  },{
    name: 'smiles'
  },{
    name: 'formula'
  },{
    name: 'predicted', type: 'bool'
  },{
    name: 'name'
  },{
    name: 'nhits', type:'number'
  },{
    name: 'atoms', defaultValue: [] // array of atom indexes of molecule which are the substructure of the query
  },{
    name: 'scans', defaultValue: [] // Filled when molecule is selected
  },{
    name: 'score', type:'number' // Only filled when scan is selected
  },{
    name: 'deltappm', type:'number' // Only filled when scan is selected
  }, {
    name: 'mim', type:'number'
  }, {
    name: 'logp', type:'number'
  }, {
    name: 'assigned', type:'bool'
  }, {
    name: 'reference'
  }]
});
