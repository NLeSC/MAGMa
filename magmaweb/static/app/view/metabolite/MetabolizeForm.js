/**
 * Windowed form to metabolize all structures.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.metabolite.MetabolizeForm', {
    extend: 'Ext.window.Window',
    title : 'Metabolize all structures',
    alias: 'widget.metabolizeform',
    requires: [
        'Ext.form.Panel',
        'Esc.magmaweb.view.metabolite.MetabolizeFieldSet',
        'Esc.magmaweb.view.fragment.AnnotateFieldSet'
    ],
    modal : true,
    height : 300,
    width : 600,
    layout : 'fit',
    closeAction : 'hide',
    items : {
        xtype : 'form',
        bodyPadding : 5,
        defaults : {
            bodyPadding : 5
        },
        border : false,
        autoScroll : true,
        items : [ {
            xtype : 'metabolizefieldset'
        }, {
            xtype : 'annotatefieldset',
            collapsed : true,
            collapsible : true
        } ],
        buttons : [ {
            text : 'Metabolize all structures',
            action: 'metabolize'
        }, {
            text : 'Reset',
            handler : function() {
                this.up('form').getForm().reset();
            }
        } ]
    },
    /**
     * Enable or disable annotate fieldset.
     *
     * @param {Boolean} disabled True to disable
     */
    setDisabledAnnotateFieldset: function(disabled) {
        this.query('annotatefieldset')[0].setDisabled(disabled);
    },
    /**
     * Get form.
     *
     * @return {Ext.form.Basic}
     */
    getForm: function() {
        return this.down('form').getForm();
    }
});