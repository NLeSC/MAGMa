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
Ext.define('Esc.magmaweb.view.molecule.AddFieldSet', {
  extend: 'Ext.form.Panel',
  alias: 'widget.addstructurefieldset',
  requires: [
    'Ext.form.field.File',
    'Esc.form.field.TextareaTab',
    'Ext.form.field.ComboBox',
    'Ext.form.field.Display',
    'Ext.tab.Panel'
  ],
  frame: true,
  bodyPadding: '5',
  items: [{
    xtype: 'tabpanel',
    activeTab: 0,
    defaults: {
      border: false,
      frame: true,
      bodyPadding: 5
    },
    plain: true,
    border: false,
    items: [{
      title: 'Database',
      id: 'structure_database_tab',
      layout: 'anchor',
      defaults: {
        labelWidth: 170
      },
      items: [{
        xtype: 'combo',
        fieldLabel: 'Chemical structure database',
        labelSeparator: '',
        name: 'structure_database',
        emptyText: 'No database selected',
        store: [
          ['pubchem', 'PubChem'],
          ['kegg', 'Kegg'],
          ['hmdb', 'HMDB']
        ],
        anchor: '80%',
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
            var excl_halo = this.up('form').down('checkbox[name="excl_halo"]');
            excl_halo.setVisible(value == 'pubchem' || value == 'kegg');
            var ids_file = this.up('form').down('filefield[name="ids_file"]');
            ids_file.setVisible(value == 'pubchem');
          }
        }
      }, {
        xtype: 'displayfield',
        value: '<br>Filter options:'
      }, {
        xtype: 'numberfield',
        fieldLabel: 'Mass',
        labelSeparator: '',
        labelAttrTpl: 'data-qwidth=200 data-qtip="Mass limit for retrieving candidate molecules"',
        afterLabelTextTpl: '<span class="relation">&le;</span>',
        name: 'max_mz',
        minValue: 1,
        value: 1200,
        anchor: '80%'
      }, {
        xtype: 'numberfield',
        fieldLabel: 'PubChem reference score',
        labelSeparator: '',
        afterLabelTextTpl: '<span class="relation">&ge;</span>',
        labelAttrTpl: 'data-qwidth=200 data-qtip="Minimum number of related PubChem substances"',
        name: 'min_refscore',
        minValue: 1,
        value: 1,
        hidden: true,
        anchor: '80%'
      }, {
        xtype: 'checkbox',
        fieldLabel: 'Excl. halogenated compounds',
        labelSeparator: '',
        labelAttrTpl: 'data-qwidth=200 data-qtip="Exclude halogenated compounds"',
        name: 'excl_halo',
        checked: true,
        hidden: true,
        anchor: '80%'
      }, {
        xtype: 'filefield',
        fieldLabel: 'Selected pubchem ID\'s',
        name: 'ids_file',
        id: 'ids_filefield',
        emptyText: 'Upload from file',
        hidden: true,
        anchor: '80%'

      }]
    }, {
      title: 'Upload',
      layout: 'anchor',
      defaults: {
        labelWidth: 150
      },
      items: [{
        xtype: 'combo',
        name: 'structure_format',
        fieldLabel: 'Format',
        store: ['smiles', 'sdf'],
        allowBlank: false,
        value: 'smiles',
        anchor: '75%'
      }, {
        xtype: 'textareatab',
        name: 'structures',
        id: 'structures_area',
        emptyText: 'Enter SDF, or smiles strings followed by space and name on each line',
        height: 200,
        anchor: '90%',
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
            molblock = molblock.replace("v2000", 'V2000');
            form.setValues({
              structure_format: 'sdf',
              structures_area: molblock
            });
          }
          return true;
        }
      }, {
        xtype: 'displayfield',
        value: 'or'
      }, {
        xtype: 'filefield',
        name: 'structures_file',
        id: 'structures_filefield',
        emptyText: 'Upload structures from file',
        anchor: '75%'
      }]
    }, {
      title: 'Draw',
      items: [{
        // during form submit fetch molblock from sketcher
        id: 'sketcher',
        contentEl: 'sketcher_content',
        border: false
      }]
    }]
  }, {
    xtype: 'metabolizefieldset',
    margin: '0 0 10 0'
  }]
});
