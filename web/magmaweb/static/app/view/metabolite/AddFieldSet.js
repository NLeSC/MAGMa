/**
 * Form element for adding structures via textarea, file upload or draw using ChemDoodle sketcher.
 *
 * The sketcher requires a dom element:
 *
 *      <div id="sketcher_content" class="x-hidden">
 *        <script type="text/javascript">
 *            var sketcher = new ChemDoodle.SketcherCanvas(
 *                'sketcher_canvas', 500, 300,
 *                '${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/icons/')}',
 *                ChemDoodle.featureDetection.supports_touch(), false);
 *            sketcher.repaint();
 *            sketcher.toolbarManager.buttonSave.disable();
 *            sketcher.toolbarManager.buttonOpen.disable();
 *        </script>
 *      </div>
 *
 *  Also ChemDoodle sketcher css and js need to be included, see http://web.chemdoodle.com/tutorial/2d-structure-canvases/sketcher-canvas
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.metabolite.AddFieldSet', {
	extend: 'Ext.form.FieldSet',
	alias : 'widget.addstructurefieldset',
	requires : ['Ext.form.field.File',
	            'Esc.form.field.TextareaTab',
	            'Ext.form.field.ComboBox',
	            'Ext.form.field.Display',
	            'Ext.tab.Panel'],
	items: [{
	    xtype: 'tabpanel',
	    activeTab : 0,
	    defaults : {
			border: false,
	        bodyPadding : 5
	    },
	    plain : true,
	    border: false,
	    items : [{
	    	title: 'Database',
	    	id: 'structure_database_tab',
		    defaults : {
		        labelWidth : 200
		    },
	    	items: [{
	    		xtype: 'combo',
                fieldLabel: '&nbsp;',
                labelSeparator: '',
	    		name: 'structure_database',
	    		emptyText: 'No database selected',
                store: [
                    ['pubchem', 'PubChem'],
                    ['kegg', 'Kegg'],
                    ['hmdb', 'Human Metabolome Database']
                ],
                listeners: {
                    /**
                     * Only show 'PubChem reference score' when 'pubchem' is selected.
                     *
                     * @param t
                     * @param value
                     */
                    change: function(t, value) {
                        var min_refscore = this.up('form').down('numberfield[name="min_refscore"]');
                        min_refscore.setVisible(value == 'pubchem');
                    }
                }
	    	}, {
    			xtype: 'numberfield',
                fieldLabel: 'PubChem reference score',
                labelSeparator: '',
                afterLabelTextTpl: '<span class="relation">&ge;</span>',
                tooltip: 'Minimum reference score based on nr. of synonyms and substances',
    			name: 'min_refscore',
    			minValue: 1,
                value: 1,
                hidden: true
    		}, {
    			xtype: 'numberfield',
                fieldLabel: 'Mass',
                labelSeparator: '',
                tooltip: 'Maximum m/z (mass-to-charge ratio)',
                afterLabelTextTpl: '<span class="relation">&le;</span>',
    			name: 'max_mz',
    			minValue: 1,
                value: 1200
	    	}]
	    }, {
	        title : 'Upload',
		    defaults : {
		        labelWidth : 200
		    },
	        items : [{
	                    xtype : 'combo',
	                    name : 'structure_format',
	                    fieldLabel: 'Format',
	                    store : ['smiles', 'sdf'],
	                    allowBlank : false,
	                    value : 'smiles'
	                }, {
	                    xtype : 'textareatab',
	                    name : 'structures',
	                    id : 'structures_area',
	                    emptyText : 'Enter smiles strings followed by space and name on each line',
	                    height : 200,
	                    width : 500,
	                    /**
	                     * Use validator to write sketched molecule in textarea as molblock.
	                     * A sketched molecule will overwrite the textarea.
	                     * @param {Object} value
	                     * @return {Boolean}
	                     */
	                    validator: function(value) {
	                        var form = this.up('form').getForm();
	                        var mol = sketcher.getMolecule();
	                        if (mol.bonds.length > 0) {
	                            var molblock = ChemDoodle.writeMOL(mol);
	                            // rdkit does not like v2000 in sdf
	                            // replace them with V2000
	                            // See https://github.com/NLeSC/MAGMa/issues/88
	                            molblock = molblock.replace("v2000",'V2000');
	                            form.setValues({
	                                structure_format: 'sdf',
	                                structures_area: molblock
	                            });
	                        }

	                        return true;
	                    }
	                }, {
	                    xtype : 'displayfield',
	                    value : 'or'
	                }, {
	                    xtype : 'filefield',
	                    name : 'structures_file',
	                    id : 'structures_filefield',
	                    emptyText : 'Upload structures from file',
	                    width : 300
	                }]
	    }, {
	        title : 'Draw',
	        items : [{
	                    id : 'sketcher', // during form submit fetch molblock
	                    // from sketcher
	                    contentEl : 'sketcher_content',
	                    border : false
	                }]
	    }]
    }]
});
