/**
 * Reaction column.
 *
 * Data example:
 *
 * {
 *   'reactants': [{
 *      'esterase': [123, 45]
 *   }]
 *   'products': [{
 *      'theogallin': [678, 90]
 *   }]
 * }
 *
 * Will be rendered as
 * <ul>
 *   <li>Reactant of<ul>
 *     <li class="reaction-col reaction-col-r-0" data-qtip="Click to see products">esterase (123/45)</li>
 *   </ul></li>
 *   <li>Product of<ul>
 *     <li class="reaction-col reaction-col-p-0" data-qtip="Click to see reactants">theogallin (678/90)</li>
 *   </ul></li>
 * </ul>
 *
 * Clicking on 'esterase' will filter molecules on which are products of 'esterase' reaction with current molecule as reactant.
 *
 */
Ext.define('Esc.magmaweb.view.metabolite.ReactionColumn', {
    extend: 'Ext.grid.column.Column',
    alias: ['widget.reactioncolumn'],
    alternateClassName: 'Ext.grid.ReactionColumn',
    defaultRenderer: function(value, meta, record) {
        // TODO convert json to html with classes
        return value;
    },
    processEvent : function(type, view, cell, recordIndex, cellIndex, e, record, row){
        var me = this,
        target = e.getTarget(),
        match,
        item, fn,
        key = type == 'keydown' && e.getKey(),
        disabled;

        // TODO attach click event to reaction and perform filtering
        // See Ext-grid-column-Action

        return me.callParent(arguments);
    }
});
