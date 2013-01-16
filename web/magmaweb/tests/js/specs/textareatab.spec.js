describe('Esc.form.field.TextareaTab', function() {

	it('tab key pressed', function() {
		var ta = Ext.create('Esc.form.field.TextareaTab', {
			renderTo: Ext.getBody()
		});

		// Mock user pressing TAB
		var event = {
			TAB: 9,
			ctrlKey: false,
			shiftKey: false,
			altKey: false,
			getKey: function() { return 9;},
			stopEvent: function() {}
		};
		ta.onSpecialKey(ta, event);

		expect(ta.getValue()).toEqual("\t");

		Ext.destroy(ta);
	});

		it('ctrl-tab key pressed', function() {
		var ta = Ext.create('Esc.form.field.TextareaTab', {
			renderTo: Ext.getBody()
		});

		// Mock user pressing ctrl-TAB
		var event = {
			TAB: 9,
			ctrlKey: true,
			shiftKey: false,
			altKey: false,
			getKey: function() { return 9;},
			stopEvent: function() {}
		};
		ta.onSpecialKey(ta, event);

		expect(ta.getValue()).toEqual("");

		Ext.destroy(ta);
	});
});