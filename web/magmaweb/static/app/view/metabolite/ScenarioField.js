/**
 * Editiable grid field metabolize scenarios.
 *
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.view.metabolite.ScenarioField', {
    extend: 'Ext.grid.Panel',
    mixins: {
        field: 'Ext.form.field.Field'
    },
    alias: 'widget.scenariofield',
    requires: [
        'Esc.magmaweb.store.Scenarios',
        'Ext.grid.column.Action',
        'Ext.selection.CellModel',
        'Ext.grid.plugin.CellEditing',
        'Ext.grid.plugin.DragDrop'
    ],
    width: 400,
    height: 140,
    initComponent: function() {
        var me = this;

        this.predefined = {
            drugs: [
                {type: 'phase1', steps: '2'},
                {type: 'phase2', steps: '1'}
            ],
            polyphenols: [
                {type: 'glycosidase', steps: 'complete'},
                {type: 'mass_filter', steps: '600'},
                {type: 'gut', steps: '6'},
                {type: 'phase2_selected', steps: '2'}
            ]
        };

        var store = Ext.create('Esc.magmaweb.store.Scenarios');

        this.editing =  Ext.create('Ext.grid.plugin.CellEditing', {
            clicksToEdit: 1
        })
        Ext.apply(this, {
           store: store,
           plugins: [this.editing],
           columns: [{
               header: 'Transformation type', dataIndex: 'type', flex:1, sortable: false, editor: {
                   xtype: 'combo',
                   typeAhead: true,
                   triggerAction: 'all',
                   store: [
                       'phase1',
                       'phase2',
                       'phase2_selected',
                       'glycosidase',
                       'mass_filter',
                       'gut'
                   ]
               }
            }, {
               header: 'Steps', dataIndex: 'steps', sortable: false, editor: {
                 xtype: 'textfield', allowBlank: false
               }
            }, {
               xtype: 'actioncolumn',
               width: 60, sortable: false,
               items: [{
                   iconCls: 'x-form-itemselector-up',
                   tooltip: 'Move transformation up',
                   scope: me,
                   handler: this.onMoveUpClick
               }, {
                   iconCls: 'x-form-itemselector-down',
                   tooltip: 'Move transformation down',
                   scope: me,
                   handler: this.onMoveDownClick
               }, {
                   iconCls: 'icon-delete',
                   tooltip: 'Delete transformation',
                   scope: me,
                   handler: this.onDeleteClick
               }]
            }],
            tbar: ['Predefined scenarios: ', {
                text: 'Drugs',
                scope: me,
                handler: me.onLoadDrugsScenario
            }, {
                text: 'Polyphenols',
                scope: me,
                handler: me.onLoadPolyphenolsScenario
            }, '-', {
                iconCls: 'icon-add',
                text: 'Add transformation',
                scope: me,
                handler: me.onAddClick
            }],
            viewConfig: {
                markDirty: false,
                emptyText: 'No scenarios defined: Add one by selecting a predefined scenario or making a custom scenario by adding transformations manually'
            },
            /**
             * @cfg {Boolean} [allowBlank=false] `false` to require at least one item in the list to be selected, `true` to allow no
             * selection.
             */
            allowBlank: false
        });

        me.callParent(arguments);
        me.initField();
    },
    setValue: function(value) {
        var me = this;
        if (value) {
            me.getStore().loadData(value);
            me.mixins.field.setValue.call(me, value);
        }
    },
    getRawValue: function() {
        var store = this.getStore();
        var data = [];
        store.data.each(function(record) {
            data.push(record.data);
        });
        return data;
    },
    /**
     * Returns the value that would be included in a standard form submit for this field.
     *
     * @return {String} The value to be submitted, or `null`.
     */
    getValue: function() {
        return Ext.JSON.encode(this.getRawValue());
    },
    onAddClick: function() {
        this.getStore().insert(0, {type: 'phase1', steps: 1});
        this.editing.startEditByPosition({row: 0, column:0});
    },
    onDeleteClick: function(grid, rowIndex) {
        this.getStore().removeAt(rowIndex);
    },
    onLoadDrugsScenario: function() {
        this.getStore().loadData(this.predefined['drugs']);
    },
    onLoadPolyphenolsScenario: function() {
        this.getStore().loadData(this.predefined['polyphenols']);
    },
    onMoveUpClick: function(grid, rowIndex) {
        var store = this.getStore();
        var index = 0;
        var record = store.getAt(rowIndex);
        index = Math.max(index, rowIndex - 1);

        store.remove(record);
        store.insert(index, record);
    },
    onMoveDownClick: function(grid, rowIndex) {
        var store = this.getStore();
        var index = store.getCount() - 1;
        var record = store.getAt(rowIndex);
        index = Math.min(index, rowIndex + 1);

        store.remove(record);
        store.insert(index, record);
    },
});