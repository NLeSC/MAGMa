/**
 * Windowed form to metabolize all structures.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.molecule.MetabolizeForm', {
    extend: 'Ext.window.Window',
    title : 'Metabolize all structures',
    alias: 'widget.metabolizeform',
    requires: [
        'Ext.form.Panel',
        'Esc.magmaweb.view.molecule.MetabolizeFieldSet',
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
        trackResetOnLoad: true,
        items : [ {
            margin: '0 0 10 0',
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
    },
    /**
     * Load form defaults from server.
     * During loading this.loading will be true.
     *
     * @param {String} url Url to load defaults from.
     */
    loadDefaults: function(url) {
        var me = this;
        this.loading = true;
        this.getForm().load({
            url: url,
            method: 'GET',
            waitMsg: 'Fetching defaults',
            failure: function(form, action) {
                delete me.loading;
                Ext.Error.raise(action.response.responseText);
            },
            success: function(form) {
                delete me.loading;
            }
        });
    }
});