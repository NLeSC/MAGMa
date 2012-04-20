/**
 * Form element for adding structures via textarea, file upload or draw using ChemDoodle sketcher.
 *
 * The sketcher requires a dom element:
 *
 *      <div id="sketcher_content" class="x-hidden">
 *        <script language="javascript">
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
    extend : 'Ext.tab.Panel',
    alias : 'widget.addstructurefieldset',
    requires : ['Ext.form.field.File', 'Ext.form.field.TextArea',
            'Ext.form.field.ComboBox', 'Ext.form.field.Display'],
    activeTab : 0,
    defaults : {
        bodyPadding : 5
    },
    plain : true,
    items : [{
        title : 'Upload',
        items : [{
                    xtype : 'combo',
                    name : 'structure_format',
                    store : ['smiles', 'sdf'],
                    allowBlank : false,
                    value : 'smiles'
                }, {
                    xtype : 'textarea',
                    name : 'structures',
                    id : 'structures_area',
                    emptyText : 'Enter smile string followed by space and name on each line',
                    height : 200,
                    width : 500,
                    /**
                     * Use validator to write sketched molecule in textarea as molblock.
                     * A sketched molecule will overwrite the textarea.
                     * @param {} value
                     * @return {Boolean}
                     */
                    validator: function(value) {
                        var form = this.up('form').getForm();
                        var mol = sketcher.getMolecule();
                        if (mol.bonds.length > 0) {
                            var molblock = ChemDoodle.writeMOL(mol);
                            form.setValues({
                                structures_format: 'sdf',
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
                    emptyText : 'Upload structures in file',
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
});