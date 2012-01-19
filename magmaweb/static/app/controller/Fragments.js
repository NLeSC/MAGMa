/**
 * Fragment controller.
 *
 * Handles actions performed on the fragment views.
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 */
Ext.define('Esc.magmaweb.controller.Fragments', {
  extend: 'Ext.app.Controller',
  stores: [ 'Fragments' ],
  models: [ 'Fragment' ],
  views: [ 'fragment.Tree' ],
  refs: [{
    ref: 'fragmentTree', selector: 'fragmenttree'
  }],
  init: function() {
    console.log('Fragments controller init');

    this.getFragmentsStore().on('load', this.onLoad, this);

    this.control({
      'fragmenttree': {
        select: this.onSelect,
        deselect: this.onDeselect,
        itemcollapse: this.onFragmentCollapse,
        itemexpand: this.onFragmentExpand
      }
    });

    this.application.on('scanandmetaboliteselect', this.loadFragments, this);
    this.application.on('scanandmetabolitenoselect', this.clearFragments, this);
    this.application.on('mspectraload', this.initMolecules, this);
    this.application.on('peakdeselect', this.clearFragmentSelection, this);
    this.application.on('peakselect', this.selectFragmentByPeak, this);

    this.application.addEvents(
      /**
       * @event
       * Triggered when a fragment node is collapsed.
       * @param {Esc.magmaweb.model.Fragment} fragment Fragment which has been collapsed.
       */
      'fragmentcollapse',
      /**
       * @event
       * Triggered when a fragment node is expanded.
       * @param {Esc.magmaweb.model.Fragment} fragment Fragment which has been expanded.
       */
      'fragmentexpand',
      /**
       * @event
       * Triggered when a fragment node is selected.
       * @param {Esc.magmaweb.model.Fragment} fragment Fragment which has been selected.
       */
      'fragmentselect',
      /**
       * @event
       * Triggered when a fragment node is deselected.
       * @param {Esc.magmaweb.model.Fragment} fragment Fragment which has been deselected.
       */
      'fragmentdeselect',
      /**
       * @event
       * Triggered when a children of a fragment node are loaded.
       * @param {Esc.magmaweb.model.Fragment} parent
       * @param {Array} children Array of fragment children.
       */
      'fragmentload'
    );
  },
  /**
   * Loads lvl 1 and 2 fragments of a metabolite scan combination.
   *
   * @param {Number} scanid Scan identifier.
   * @param {Number} metid Metabolite idenfitier.
   */
  loadFragments: function (scanid, metid) {
    this.clearFragments();
    console.log('Show fragments of scan '+scanid+' metabolite '+metid);
    var store = this.getFragmentsStore();
    store.setProxy(this.fragmentProxyFactory(scanid, metid));
    store.load();
  },
  /**
   * Need to change url of fragment proxy so use a factory to create a new proxy foreach scan/metabolite combo
   *
   * @param {Number} scanid Scan identifier.
   * @param {Number} metid Metabolite idenfitier.
   * @private
   */
  fragmentProxyFactory: function (scanid, metid) {
    return Ext.create('Ext.data.proxy.Ajax', {
      // url is build when scan and metabolite are selected
      url: Ext.String.format(this.application.getUrls().fragments, scanid, metid),
      reader: {
          type: 'json',
          root: 'children',
          idProperty: 'fragid'
      }
    });
  },
  /**
   * Clears fragments from store.
   */
  clearFragments: function() {
    console.log('Clearing fragments and mspectra >lvl1');
    this.getFragmentsStore().getRootNode().removeAll();
  },
  onFragmentCollapse: function(fragment) {
    this.application.fireEvent('fragmentcollapse', fragment);
  },
  onFragmentExpand: function(fragment) {
    if (fragment.firstChild == null) {
      return; // root node auto expands, but is no fragment, so dont fire event
    }
    this.application.fireEvent('fragmentexpand', fragment);
  },
  onSelect: function(rm, r) {
    console.log('Selected fragment '+r.id);
    // show child mspectra of selected node or mz
    if (!r.isLeaf()) {
      // onselect then expand
      if (r.isExpanded()) {
        this.onFragmentExpand(r);
      } else {
        r.expand();
      }
    }
    this.application.fireEvent('fragmentselect', r);
  },
  onDeselect: function(rm, fragment) {
      this.application.fireEvent('fragmentdeselect', fragment);
  },
  /**
   * Clears fragment selection.
   */
  clearFragmentSelection: function() {
      this.getFragmentTree().getSelectionModel().deselectAll();
  },
  onLoad: function(t, parent, children) {
    this.application.fireEvent('fragmentload', parent, children);
  },
  selectFragment: function(fragment) {
    this.getFragmentTree().getSelectionModel().select([fragment]);
    if (!fragment.isLeaf()) {
      if (fragment.isExpanded()) {
        this.application.fireEvent('fragmentexpand', fragment);
      } else {
        fragment.expand();
      }
    }
  },
  /**
   * When user selects peak in spectra then select the fragment belonging to peak in fragment tree
   *
   * @param {Number} mz m/z of peak
   * @param {Number} mslevel MS level of peak
   */
  selectFragmentByPeak: function(mz, mslevel) {
    // find fragment based on mz + mslevel
    var node = this.getFragmentsStore().getNodeByMzMslevel(mz, mslevel);
    this.selectFragment(node);
  },
  /**
   * Forces molecules canvases to be drawn
   */
  initMolecules: function() {
    this.getFragmentTree().initMolecules();
  }
});
