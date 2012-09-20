/**
 * Form to upload ms data.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.scan.UploadForm', {
    extend: 'Ext.form.Panel',
    alias: 'widget.scanuploadform',
    autoScroll: true,
    trackResetOnLoad: true,
    items: [{
        xtype: 'uploadmsdatafieldset'
    }, {
        xtype : 'annotatefieldset',
        collapsed : true,
        collapsible : true
    }],
    buttons: [{
        text: 'Upload ms data',
        action: 'uploadmsdata'
    }, {
        text: 'Reset',
        handler: function() {
            this.up('form').getForm().reset();
        }
    }, {
       text: 'Cancel',
       action: 'uploadmsdatacancel'
    }],
    /**
     * Enable or disable annotate fieldset.
     *
     * @param {Boolean} disabled True to disable
     */
    setDisabledAnnotateFieldset: function(disabled) {
        this.query('annotatefieldset')[0].setDisabled(disabled);
    },
    /**
     * Load form defaults from server.
     *
     * @param {String} url Url to load defaults from.
     */
    loadDefaults: function(url) {
        this.load({
            url: url,
            method: 'GET',
            waitMsg: 'Fetching defaults',
            failure: function(form, action) {
                Ext.Error.raise(action.response.responseText);
            }
        });
    }
});