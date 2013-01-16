/**
 * TextArea which enteres TAB key instead of changing focus to next form field.
 */
Ext.define('Esc.form.field.TextareaTab', {
	extend: 'Ext.form.field.TextArea',
	alias : 'widget.textareatab',
	initComponent: function() {
		this.callParent(arguments);
		this.on('specialkey', this.onSpecialKey, this);
	},
	/**
	 * @private Event handle fired when the user presses special key like TAB.
	 */
	onSpecialKey: function(field, event) {
		if (event.getKey() == event.TAB
		&& !event.shiftKey
		&& !event.altKey
		&& !event.ctrlKey
		) {
			var i = field.inputEl.dom;
			var v = i.value;
			var start = i.selectionStart;
			var end = i.selectionEnd;
			i.value = v.substring(0, start) + "\t" + v.substring(end);
			i.selectionStart = i.selectionEnd = start + 1;

			event.stopEvent();
		}
	}
});