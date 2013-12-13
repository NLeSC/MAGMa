/**
 * Reaction column.
 *
 * Data example:
 *
 * {
 *   'reactantof': [{
 *      'esterase': {'nr': 123, 'nrp': 45}
 *   }]
 *   'productof': [{
 *      'theogallin': {'nr: 678, 'nrp': 90}
 *   }]
 * }
 *
 * Will be rendered as
 * <ul>
 *   <li class="reaction-col reaction-col-r>Reactant of<ul>
 *     <li class="reaction-col reaction-col-r-0" data-qtip="Click to see products">esterase (123/45)</li>
 *   </ul></li>
 *   <li class="reaction-col reaction-col-p>Product of<ul>
 *     <li class="reaction-col reaction-col-p-0" data-qtip="Click to see reactants">theogallin (678/90)</li>
 *   </ul></li>
 * </ul>
 *
 * Clicking on 'esterase' will filter molecules on which are products of 'esterase' reaction with current molecule as reactant.
 *
 */
Ext.define('Esc.magmaweb.view.metabolite.ReactionColumn', {
    extend: 'Ext.grid.column.Template',
    alias: ['widget.reactioncolumn'],
    alternateClassName: 'Ext.grid.ReactionColumn',
    tpl: '<ul>' +
         '  <li class="reaction-col-r-g"><span class="reaction-col-item" data-qtip="Click to see products of this molecule">Reactant of</span><ul class="reaction-col-items">' +
         '    <tpl foreach="reactionsequence.reactantof"><li>- <span class="reaction-col-item reaction-col-r-{[xindex]} data-qtip="Click to see products of this reaction">{$}</span> ({nr}/{nrp})</li></tpl>' +
         '  </ul></li>' +
         '  <li class="reaction-col-p-g"><span class="reaction-col-item" data-qtip="Click to see reactants of this molecule">Product of</span><ul class="reaction-col-items">' +
         '    <tpl foreach="reactionsequence.productof"><li>- <span class="reaction-col-item reaction-col-p-{[xindex]}" data-qtip="Click to see reactants of this reaction">{$}</span> ({nr}/{nrp})</li></tpl>' +
         '  </ul></li>' +
         '</ul>',
    groupRe: new RegExp('reaction-col-(\\w)-g'),
    itemRe: new RegExp('reaction-col-(\\w)-(\\d+)'),
    processEvent : function(type, view, cell, recordIndex, cellIndex, e, record, row){
        var me = this,
        target = e.getTarget(),
        match,
        item, fn,
        key = type == 'keydown' && e.getKey(),
        disabled;
        var reaction_filter = view.up('grid').filters.getFilter('reactionsequence');

//        if (key && !Ext.fly(target).findParent(view.getCellSelector())) {
//            target = Ext.fly(cell).down('.reaction-col-item', true);
//        }
//
//        if (type === 'click') {
//            // reactant of or product of clicked
//            match = target.className.match(me.groupRe);
//            if match {
//                if (match[1] === 'r') {
//                    reaction_filter.setValue({
//                        'reactant': record.data.metid
//                    });
//                } else if (match[1] === 'p') {
//                    reaction_filter.setValue({
//                        'product': record.data.metid
//                    });
//                }
//            }
//
//            // reaction name
//            match = target.className.match(me.itemRe);
//            if match {
//                if (match[1] === 'r') {
//                    reaction_filter.setValue({
//                        'reactant': record.data.metid,
//                    });
//                } else if (match[1] === 'p') {
//                    reaction_filter.setValue({
//                        'product': record.data.metid
//                    });
//                }
//            }
//        }
        // TODO attach click event to reaction and perform filtering
        // See Ext-grid-column-Action

        return me.callParent(arguments);
    }
});
