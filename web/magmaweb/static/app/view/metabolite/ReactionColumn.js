/**
 * Reaction column.
 *
 * Data example:
 *
 * {
 *   'reactantof': [{
 *      'esterase': {'nr': 123, 'nrp': 45}
 *   }],
 *   'productof': [{
 *      'theogallin': {'nr': 678, 'nrp': 90}
 *   }]
 * }
 *
 * Will be rendered as
 * <ul>
 *     <li><span class="reaction-col-item reaction-col-r-g" data-qtip="Click to see products of this molecule"><U>Products via</U></span>
 *         <ul class="reaction-col-items">
 *            <li>&nbsp&nbsp&nbsp&nbsp<span class="reaction-col-item reaction-col-r-i" data-qtip="Click to see products of this molecule via esterase reaction">esterase</span> (123/45)</li>
 *         </ul>
 *     </li>
 *     <li><span class="reaction-col-item reaction-col-p-g" data-qtip="Click to see reactants of this molecule"><U>Reactants via</U></span>
 *        <ul class="reaction-col-items">
 *            <li>&nbsp&nbsp&nbsp&nbsp<span class="reaction-col-item reaction-col-p-i" data-qtip="Click to see reactants of this molecule via theogallin reaction">theogallin</span> (678/90)</li>
 *        </ul>
 *    </li>
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
         '  <tpl if="reactionsequence.reactants">' +
         '    <li><U><span class="reaction-col-item reaction-col-p-g" data-qwidth="100" data-qtip="Click to see reactants of this molecule">Reactants via</span></U><ul class="reaction-col-items">' +
         '    <table border=0>' +
         '    <tpl foreach="reactionsequence.reactants">'+
         '      <tr><td style="width:18px">' +
         '      <tpl if="nrp"> <i><b>{nr}</i></b>' +
         '      <tpl else>{nr}' +
         '      </tpl></td>' +
         '      <td><span class="reaction-col-item reaction-col-p-i" data-qwidth="100" data-qtip="Click to see reactants of this molecule via {$} reaction">{$}</span></td></tr>' +
         '    </tpl></table>' +
         '  </ul></li>' +
         '  </tpl>' +
         '  <tpl if="reactionsequence.products">' +
         '    <p style="margin:3px 0;"></p>' +
         '    <li><U><span class="reaction-col-item reaction-col-r-g" data-qwidth="100" data-qtip="Click to see products of this molecule">Products via</span></U><ul class="reaction-col-items">' +
         '    <table border=0>' +
         '    <tpl foreach="reactionsequence.products">' +
         '      <tr><td style="width:18px">' +
         '      <tpl if="nrp"> <i><b>{nr}</i></b>' +
         '      <tpl else>{nr}' +
         '      </tpl></td>' +
         '      <td><span class="reaction-col-item reaction-col-r-i" data-qwidth="100" data-qtip="Click to see products of this molecule via {$} reaction">{$}</span></td></tr>' +
         '    </tpl></table>' +
         '  </ul></li>' +
         '  </tpl>' +
         '</ul>',
    itemRe: new RegExp('reaction-col-(\\w)-(\\w)'),
    getFilter: function(view) {
        return view.up('grid').filters.getFilter('reactionsequence');
    },
    processEvent : function(type, view, cell, recordIndex, cellIndex, e, record, row){
        var me = this, target = e.getTarget(), match;
        var reaction_filter = this.getFilter(view);

        // only listen for click events
        if (type !== 'click') {
            return false;
        }

        if (match = target.className.match(me.itemRe)) {
            if (match[2] === 'g') {
                if (match[1] === 'r') {
                    // 'reactant of' clicked
                    reaction_filter.setValue({
                        'reactant': record.data.metid
                    });
                    return reaction_filter.setActive(true);
                } else if (match[1] === 'p') {
                    // 'product of' clicked
                    reaction_filter.setValue({
                        'product': record.data.metid
                    });
                    return reaction_filter.setActive(true);
                }
            } else if (match[2] === 'i') {
                if (match[1] === 'r') {
                    // reaction name of reactant clicked
                    reaction_filter.setValue({
                        'name': target.innerText,
                        'reactant': record.data.metid,
                    });
                    return reaction_filter.setActive(true);
                } else if (match[1] === 'p') {
                    // reaction name of product clicked
                    reaction_filter.setValue({
                        'name': target.innerText,
                        'product': record.data.metid
                    });
                    return reaction_filter.setActive(true);
                }
            }
        }

        // do not select molecule with this column.
        return false;
    }
});
