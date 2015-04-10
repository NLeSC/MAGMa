/**
 * Reaction filter with 3 fields:
 * <ol>
 * <li>Reactant, molecule identifier of reactant as numeric field</li>
 * <li>Product, molecule identifier of product as numeric field</li>
 * <li>Name, Name of reaction as a string field</li>
 * </ol>
 *
 * Usages:
 * <ul>
 * <li>Show molecules which have current molecule as a reactant</li>
 * <li>Show molecules which have current molecule as a product</li>
 * <li>Show molecules which are reactants of a reaction (eg. esterase)</li>
 * <li>Show molecules which are products of a reaction (eg. esterase)</li>
 * <li>Show molecules which are reactants of a esterase reaction with current molecule as product</li>
 * <li>Show molecules which are products of a esterase reaction with current molecule as reactant</li>
 * </ul>
 * Current molecule is the molecule filtered on as reactant or product.
 *
 */
Ext.define('Esc.magmaweb.view.molecule.ReactionFilter', {
  extend: 'Ext.ux.grid.filter.Filter',
  alias: ['gridfilter.reaction'],
  createMenu: function(config) {
    var me = this;
    var menu = Ext.create('Ext.menu.Menu');
    me.reactant = Ext.create('Ext.form.field.Number', {
      enableKeyEvents: true,
      labelWidth: 120,
      fieldLabel: 'Products of (ID)',
    });
    me.reactant.on('specialkey', this.onSpecialKey, this);
    menu.add(me.reactant);

    me.product = Ext.create('Ext.form.field.Number', {
      enableKeyEvents: true,
      labelWidth: 120,
      fieldLabel: 'Reactants of (ID)'
    });
    me.product.on('specialkey', this.onSpecialKey, this);
    menu.add(me.product);

    me.reaction_name = Ext.create('Ext.form.field.Text', {
      enableKeyEvents: true,
      labelWidth: 120,
      fieldLabel: 'Via reaction'
    });
    me.reaction_name.on('specialkey', this.onSpecialKey, this);
    menu.add(me.reaction_name);

    return menu;
  },
  getValue: function() {
    var me = this;
    var value = {};
    var v;
    if (me.reactant.isValid()) {
      v = me.reactant.getValue();
      if (v !== null && v !== "") {
        value.reactant = v;
      }
    }
    if (me.product.isValid()) {
      v = me.product.getValue();
      if (v !== null && v !== "") {
        value.product = v;
      }
    }
    if (me.reaction_name.isValid()) {
      v = me.reaction_name.getValue();
      if (v !== null && v !== "") {
        value.name = v;
      }
    }
    return value;
  },
  setValue: function(value) {
    var me = this;
    if ('reactant' in value) {
      me.reactant.setValue(value.reactant);
    } else {
      me.reactant.reset();
    }
    if ('product' in value) {
      me.product.setValue(value.product);
    } else {
      me.product.reset();
    }
    if ('name' in value) {
      me.reaction_name.setValue(value.name);
    } else {
      me.reaction_name.reset();
    }
    this.fireUpdate();
  },
  getSerialArgs: function() {
    var args = this.getValue();
    args.type = 'reaction';
    return args;
  },
  validateRecord: function(record) {
    // TODO implement
    return true;
  },
  onSpecialKey: function(field, e) {
    if (e.getKey() === e.ENTER) {
      this.fireUpdate();
    }
  }
});
