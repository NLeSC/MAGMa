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
  uses: [ 'Esc.magmaweb.view.fragment.AnnotateForm' ],
  /**
   * Can only annotate when there are structures and ms data.
   * @property {Object} annotatable
   * @property {Boolean} annotabable.structures Whether there are structures
   * @property {Boolean} annotabable.msdata Whether there is ms data
   */
  annotatable: { structures: false, msdata: false },
  init: function() {
    console.log('Fragments controller init');

    this.getFragmentsStore().on('load', this.onLoad, this);

    this.control({
      'fragmenttree': {
        select: this.onSelect,
        deselect: this.onDeselect,
        itemcollapse: this.onFragmentCollapse,
        itemexpand: this.onFragmentExpand
      },
      'annotateform component[action=annotate]': {
          click: this.annotateHandler
      },
      '#annotateaction': {
          click: this.annotateAction
      }
    });

    this.application.on('scanandmetaboliteselect', this.loadFragments, this);
    this.application.on('scanandmetabolitenoselect', this.clearFragments, this);
    this.application.on('mspectraload', this.initMolecules, this);
    this.application.on('peakdeselect', this.clearFragmentSelection, this);
    this.application.on('peakselect', this.selectFragmentByPeak, this);

    this.application.on('metaboliteload', function(store) {
        this.annotatable.structures = store.getTotalUnfilteredCount() > 0;
        if (this.annotatable.structures && this.annotatable.msdata) {
            this.getAnnotateActionButton().enable();
        } else {
            this.getAnnotateActionButton().disable();
        }
    }, this);
    this.application.on('chromatogramload', function(chromatogram) {
        this.annotatable.msdata = chromatogram.data.length > 0;
        if (this.annotatable.structures && this.annotatable.msdata) {
            this.getAnnotateActionButton().enable();
        } else {
            this.getAnnotateActionButton().disable();
        }
    }, this);
    this.application.on('rpcsubmitsuccess', this.rpcSubmitted, this);

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
      listeners: {
        exception: function(proxy, response, operation) {
          Ext.Error.raise({
            msg: 'Failed to load fragments from server',
            response: response,
            operation: operation
          });
        }
      },
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
    if (fragment.firstChild === null) {
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
    // when parent is root node then expand it
    if (parent.isRoot()) {
        parent.expand();
    }
    this.application.fireEvent('fragmentload', parent, parent.childNodes);
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
  },
  /**
   * Shows annotate form in modal window
   */
  showAnnotateForm: function() {
    var me = this;
    if (!this.annotateForm) {
        this.annotateForm = Ext.create('Esc.magmaweb.view.fragment.AnnotateForm');
        this.annotateForm.loadDefaults(me.application.runInfoUrl());
    }
    this.annotateForm.show();
  },
  /**
   * Handler for submission of annotate form.
   */
  annotateHandler: function() {
    var me = this;
    var wf = this.annotateForm;
    var form = wf.getForm();
    if (form.isValid()) {
      form.submit({
        url: this.application.rpcUrl('annotate'),
        waitMsg: 'Submitting action ...',
        submitEmptyText: false,
        success: function(fp, o) {
          var response = Ext.JSON.decode(o.response.responseText);
          me.application.fireEvent('rpcsubmitsuccess', response.jobid);
          wf.hide();
        },
        failure: function(form, action) {
            wf.hide();
            if (action.failureType === "server") {
              Ext.Error.raise(Ext.JSON.decode(action.response.responseText));
            } else {
              Ext.Error.raise(action.response.responseText);
            }
        }
      });
    }
  },
  getAnnotateActionButton: function() {
      return Ext.getCmp('annotateaction');
  },
  /**
   * Pols status of jobid and when completed redirects to results of new job.
   * @param {String} jobid Identifier of new job
   */
  rpcSubmitted: function(jobid) {
      var me = this;
      var app = this.application;
      /**
       * @property {String} newjobid
       * Keep track of id of submitted job
       */
      me.newjobid = jobid;
      // Overwrite annotate button to waiting/cancel button
      var annot_button = this.getAnnotateActionButton();
      annot_button.setIconCls('icon-loading');
      annot_button.setText('Waiting');
      annot_button.setTooltip('Job submitted, waiting for completion');
      annot_button.setHandler(function() {
          Ext.MessageBox.confirm('Cancel job', 'Job is still running. Do you want to cancel it?', function(but) {
              if (but == 'yes') {
                  // TODO cancel job
                  console.log('Cancelling job');
              }
          });
      });
      annot_button.enable();
      /**
       * @property {Ext.util.TaskRunner.Task} pollTask
       * Polls status of submitted job
       */
      me.pollTask = Ext.TaskManager.start({
        run: this.pollJobStatus,
        interval: 5000,
        scope: this
      });
  },
  /**
   * Polls status of job with id this.newjobid on server.
   *
   * When status is STOPPED then change annotate button to fetch result button.
   */
  pollJobStatus: function() {
    var me = this;
    var app = this.application;
    Ext.Ajax.request({
      url: app.urls.home+'status/'+me.newjobid+'.json',
      success: function(o) {
        var response = Ext.JSON.decode(o.responseText);
        console.log(response.status);
        if (response.status == 'STOPPED') {
          Ext.TaskManager.stop(me.pollTask);
          delete me.pollTask;
          var annot_button = me.getAnnotateActionButton();
          annot_button.setIconCls('');
          annot_button.setText('Fetch result');
          annot_button.setTooltip('Job completed, fetch results');
          annot_button.setHandler(function() {
              Ext.MessageBox.confirm('Fetch result', 'Job has been completed. Do you want to fetch results?', function(but) {
                  if (but == 'yes') {
                      window.location = app.urls.home+'results/'+me.newjobid;
                  }
              });
          });
          annot_button.enable();
        } else {
            me.getAnnotateActionButton().setTooltip('Job '+response.status+', waiting for completion');
        }
      },
      failure: function() {
        Ext.TaskManager.stop(me.pollTask);
        delete me.pollTask;
      },
      scope: me
    });
  },
  annotateAction: function() {
      if ('newjobid' in this) {
          // job running or completed, do not annotate
      } else {
          this.showAnnotateForm();
      }
  }
});
