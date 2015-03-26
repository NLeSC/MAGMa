describe('Fragments', function() {
  describe('store', function() {
    var url = 'data/fragments.s1133.m352.json';
    var store = null;

    beforeEach(function() {
      if (!store) {
        store = Ext.create('Esc.magmaweb.store.Fragments');
      }
    });

    it('create', function() {
      expect(store).toBeDefined();
    });

    it('getById', function() {
      var fragid = 2469;
      spyOn(store, 'getNodeById');

      store.getById(fragid);
      expect(store.getNodeById).toHaveBeenCalledWith(fragid);
    });

    it('getNodeByMzMslevel', function() {
      store.setProxy(
        Ext.create('Ext.data.proxy.Ajax', {
          // url is build when scan and molecule are selected
          url: url,
          reader: {
              type: 'json',
              root: 'children',
              idProperty: 'fragid'
          }
        })
      );
      store.load();
      waitsFor(
        function() { return !store.isLoading();},
        'Fragment store never loaded',
        5000
      );

      runs(function() {
        var node = store.getNodeByMzMslevel(122.0373001, 2);
        expect(node.getId()).toEqual(2471);
      });
    });
  });

  describe('controller', function() {
    var store = null, ctrl= null;
    var assignbut = null;

    beforeEach(function() {
       if (!ctrl) {
           ctrl = Application.getController('Fragments');
       }

       if (!store) {
           store = ctrl.getStore('Fragments');
       }
    });

    /**
     * Loads fragment tree with static json file
     */
    var fill = function(callback) {
      assignbut = jasmine.createSpyObj('abut', [ 'setParams', 'disable', 'enable', 'toggle']);
      spyOn(ctrl,'getAssignStruct2PeakButton').andReturn(assignbut);
      ctrl.loadFragments(1133, 352);
      waitsFor(
        function() { return !store.isLoading();},
        'Fragment store never loaded',
        5000
      );
      runs(callback);
    }

    it('loadFragments', function() {
      fill(function() {
        expect(store.getRootNode().hasChildNodes()).toBeTruthy();
        expect(store.getById(2469)).toBeDefined();
        expect(store.getById(2471)).toBeDefined();
      });
    });

    it('clearFragments', function() {
      // first fill then clear
      fill(function() {
        expect(store.getRootNode().hasChildNodes()).toBeTruthy();
        ctrl.clearFragments();
        expect(store.getRootNode().hasChildNodes()).toBeFalsy();
        expect(assignbut.disable).toHaveBeenCalled();
        expect(assignbut.toggle).toHaveBeenCalledWith(false);
      });
    });

    it('onFragmentExpand', function() {
      fill(function() {
        var f = { callback: function() {} };
        spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
        Ext.util.Observable.capture(ctrl.application, f.callback);

        var frag = store.getById(2469);
        frag.collapse();
        expect(frag.isExpanded()).toBeFalsy();

        frag.expand();
        ctrl.onFragmentExpand(frag);

        expect(f.callback).toHaveBeenCalledWith('fragmentexpand', jasmine.any(Object));
        var expandedfrag = f.callback.mostRecentCall.args[1];
        expect(expandedfrag.isExpanded()).toBeTruthy();
        Ext.util.Observable.releaseCapture(ctrl.application);
      });
    });

    it('onFragmentCollapse', function() {
      fill(function() {
        var f = { callback: function() {} };
        spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
        Ext.util.Observable.capture(ctrl.application, f.callback);

        var frag = store.getById(2469);
        expect(frag.isExpanded()).toBeTruthy();

        frag.collapse();
        ctrl.onFragmentCollapse(frag);
        expect(frag.isExpanded()).toBeFalsy();

        expect(f.callback).toHaveBeenCalledWith('fragmentcollapse', jasmine.any(Object));
        var collapsedfrag = f.callback.mostRecentCall.args[1];
        expect(collapsedfrag.isExpanded()).toBeFalsy();
        Ext.util.Observable.releaseCapture(ctrl.application);
      });
    });

    describe('onSelect', function() {
      it('on leaf fragment', function() {
        fill(function() {
          var frag = store.getById(2470);
          var rm = Ext.create('Ext.selection.RowModel');

          var f = { callback: function() {} };
          spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
          Ext.util.Observable.capture(ctrl.application, f.callback);

          ctrl.onSelect(rm, frag);

          expect(f.callback).toHaveBeenCalledWith('fragmentselect', jasmine.any(Object));
          Ext.util.Observable.releaseCapture(ctrl.application);
        });
      });

      it('on expanded fragment', function() {
        fill(function() {
          var frag = store.getById(2471);
          var rm = Ext.create('Ext.selection.RowModel');

          var f = { callback: function() {} };
          spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
          Ext.util.Observable.capture(ctrl.application, f.callback);

          ctrl.onSelect(rm, frag);

          expect(f.callback).toHaveBeenCalledWith('fragmentexpand', jasmine.any(Object));
          expect(f.callback).toHaveBeenCalledWith('fragmentselect', jasmine.any(Object));
          Ext.util.Observable.releaseCapture(ctrl.application);
        });
      });

      it('on collapsed fragment', function() {
        fill(function() {
          var frag = store.getById(2471);
          spyOn(frag, 'expand');
          var rm = Ext.create('Ext.selection.RowModel');

          var f = { callback: function() {} };
          spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
          Ext.util.Observable.capture(ctrl.application, f.callback);

          ctrl.onSelect(rm, frag);

          expect(frag.expand).toHaveBeenCalled();
          expect(f.callback).toHaveBeenCalledWith('fragmentselect', jasmine.any(Object));
          Ext.util.Observable.releaseCapture(ctrl.application);
        });
      });

      it('should not select lvl1 fragment', function() {
        fill(function() {
          var frag = store.getById(2469);
          var rm = Ext.create('Ext.selection.RowModel');

          var f = { callback: function() {} };
          spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
          Ext.util.Observable.capture(ctrl.application, f.callback);

          ctrl.onSelect(rm, frag);

          expect(f.callback).not.toHaveBeenCalled();
          Ext.util.Observable.releaseCapture(ctrl.application);
        });
      });
    });

    describe('onDeselect', function() {
      it('should fire event', function() {
	      var f = { callback: function() {} };
	      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
	      Ext.util.Observable.capture(ctrl.application, f.callback);

	      var frag = store.getById(2470);
	      var rm = Ext.create('Ext.selection.RowModel');

	      ctrl.onDeselect(rm, frag);

	      expect(f.callback).toHaveBeenCalledWith('fragmentdeselect', jasmine.any(Object));
	      Ext.util.Observable.releaseCapture(ctrl.application);
      });

      it('should not fire event when lvl1 fragment is unselected', function() {
	      var f = { callback: function() {} };
	      spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
	      Ext.util.Observable.capture(ctrl.application, f.callback);

	      var frag = store.getById(2469);
	      var rm = Ext.create('Ext.selection.RowModel');

	      ctrl.onDeselect(rm, frag);

	      expect(f.callback).not.toHaveBeenCalled();
	      Ext.util.Observable.releaseCapture(ctrl.application);
      });
    });

    it('clearFragmentSelection', function() {
      // mock tree
      var sm = { deselectAll: function() {} };
      var tree = { getSelectionModel: function() { return sm; } };
      spyOn(ctrl, 'getFragmentTree').andReturn(tree);
      spyOn(sm, 'deselectAll');

      ctrl.clearFragmentSelection();

      expect(sm.deselectAll).toHaveBeenCalled();
    });

    describe('onLoad', function() {
      it('root node', function() {
        var data = {
          molid: 123,
          scanid: 45,
          mz: 6789,
          isAssigned: true
        };
        var node = {
          isRoot: function() { return true;},
          expand: function() {},
          childNodes: [{
            data: data
          }]
        };
        spyOn(node, 'expand');
        var f = { callback: function() {} };
        spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
        Ext.util.Observable.capture(ctrl.application, f.callback);
        assignbut = jasmine.createSpyObj('abut', [ 'setParams', 'disable', 'enable', 'toggle']);
        spyOn(ctrl, 'getAssignStruct2PeakButton').andReturn(assignbut);

        store.fireEvent('load', store, node, 'bar');

        expect(assignbut.setParams).toHaveBeenCalledWith({
          molid: 123,
          scanid: 45,
          mz: 6789});
        expect(assignbut.enable).toHaveBeenCalled();
        expect(assignbut.toggle).toHaveBeenCalledWith(true);
        expect(f.callback).toHaveBeenCalledWith('fragmentload', node, [{data: data}]);
        expect(node.expand).toHaveBeenCalled();
        Ext.util.Observable.releaseCapture(ctrl.application);
      });

      it('non root node', function() {
          var node = {
            isRoot: function() { return false;},
            expand: function() {},
            childNodes: 'foo'
          };
          spyOn(node, 'expand');
          var f = { callback: function() {} };
          spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
          Ext.util.Observable.capture(ctrl.application, f.callback);

          store.fireEvent('load', store, node, 'bar');

          expect(f.callback).toHaveBeenCalledWith('fragmentload', node, 'foo');
          expect(node.expand).not.toHaveBeenCalled();
          Ext.util.Observable.releaseCapture(ctrl.application);
      });
    });

    describe('selectFragment', function() {
      beforeEach(function() {
        // mock tree
        var sm = { select: function() {} };
        var tree = {
        	getSelectionModel: function() { return sm; },
        	setLoading: function() {}
        };
        spyOn(ctrl, 'getFragmentTree').andReturn(tree);
        spyOn(sm, 'select');
      })

      it('leaf', function() {
        fill(function() {
          var f = { callback: function() {} };
          spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
          Ext.util.Observable.capture(ctrl.application, f.callback);

          var frag = store.getById(2470);
          ctrl.selectFragment(frag);

          expect(f.callback).not.toHaveBeenCalled();
          Ext.util.Observable.releaseCapture(ctrl.application);
        });
      });

      it('expanded', function() {
        fill(function() {
          var f = { callback: function() {} };
          spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
          Ext.util.Observable.capture(ctrl.application, f.callback);

          var frag = store.getById(2469);
          ctrl.selectFragment(frag);

          expect(f.callback).toHaveBeenCalledWith('fragmentexpand', jasmine.any(Object));
          Ext.util.Observable.releaseCapture(ctrl.application);
        });
      });

      it('collapsed', function() {
        fill(function() {
          var f = { callback: function() {} };
          spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
          Ext.util.Observable.capture(ctrl.application, f.callback);

          var frag = store.getById(2469);
          spyOn(frag, 'expand');
          frag.collapse();
          ctrl.selectFragment(frag);

          expect(frag.expand).toHaveBeenCalled();
          Ext.util.Observable.releaseCapture(ctrl.application);
        });
      });
    });

    describe('selectFragmentByPeak', function() {
    	it('should select fragment when peak has fragment', function() {
	      fill(function() {
	          spyOn(ctrl, 'selectFragment');
	          ctrl.selectFragmentByPeak(122.0373001, 2);
	          expect(ctrl.selectFragment).toHaveBeenCalled();
	          expect(ctrl.selectFragment.mostRecentCall.args[0].data.fragid).toEqual(2471);
          });
    	});

    	it('should not select fragment when peak has no fragment does not exist', function() {
	      fill(function() {

	    	  spyOn(ctrl, 'selectFragment');

	          ctrl.selectFragmentByPeak(1193.35278320312, 1);

	          expect(ctrl.selectFragment).not.toHaveBeenCalled();
          });
    	});
    });

    it('initMolecules', function() {
      // mock tree
      var tree = { initMolecules: function() {} };
      spyOn(ctrl, 'getFragmentTree').andReturn(tree);
      spyOn(tree, 'initMolecules');

      ctrl.initMolecules();

      expect(tree.initMolecules).toHaveBeenCalled();
    });

    it('proxy exception', function() {
        spyOn(Ext.Error, 'handle').andReturn(true);

        var proxy = ctrl.fragmentProxyFactory(2,3);
        proxy.fireEvent('exception', proxy, 'bla', 'foo');

        expect(Ext.Error.handle).toHaveBeenCalledWith({
            msg: 'Failed to load fragments from server',
            response: 'bla',
            operation: 'foo'
        });
    });

    it('showAnnotateForm', function() {
        ctrl.showAnnotateForm();

        waitsFor(
            function() { return !ctrl.annotateForm.loading;},
            'Form defaults never loaded',
            1000
        );

        runs(function() {
            expect(ctrl.annotateForm.getForm().getValues()).toEqual({
                ionisation_mode: -1,
                max_broken_bonds: '3',
                max_water_losses: '1',
                precursor_mz_precision: '0.001',
                ms_intensity_cutoff: '200000',
                msms_intensity_cutoff: '10',
                mz_precision: '5',
                mz_precision_abs: '0.001'
            });
        });
        expect(ctrl.annotateForm.isVisible()).toBeTruthy();
        ctrl.annotateForm.hide();
    });

    it('annotateHandler', function() {
        ctrl.showAnnotateForm();
        var form = ctrl.annotateForm.getForm();
        spyOn(form, 'submit');

        ctrl.annotateHandler();

        expect(form.submit).toHaveBeenCalledWith({
            url: '/rpc/'+Application.jobid+'/annotate',
            submitEmptyText : false,
            waitMsg: jasmine.any(String),
            success: jasmine.any(Function),
            failure: jasmine.any(Function)
        });
        ctrl.annotateForm.destroy();
    });

    it('getAnnotateActionButton', function() {
      spyOn(Ext,'getCmp');
      ctrl.getAnnotateActionButton();
      expect(Ext.getCmp).toHaveBeenCalledWith('annotateaction');
    });

    it('rpcSubmitted', function() {
      var newjobid = '3ad25048-26f6-11e1-851e-00012e260791';
      spyOn(Ext.TaskManager, 'start');
      button = {
        setIconCls: function() {},
        setText: function() {},
        setTooltip: function() {},
        setHandler: function() {},
        enable: function() {}
      };
      spyOn(ctrl, 'getAnnotateActionButton').andReturn(button);

      ctrl.rpcSubmitted(newjobid);

      expect(ctrl.newjobid).toEqual(newjobid);
      expect(Ext.TaskManager.start).toHaveBeenCalledWith({
          run: ctrl.pollJobStatus,
          interval: 5000,
          scope: ctrl
      });
    });

    describe('pollJobStatus', function() {
        var annotateActionButton;

        beforeEach(function() {
            annotateActionButton = {
                setTooltip: function() {},
                setIconCls: function() {},
                setText: function() {},
                setHandler: function() {},
                enable: function() {}
            };
        });

        it('request', function() {
            ctrl.newjobid = '3ad25048-26f6-11e1-851e-00012e260791';
            spyOn(Ext.Ajax, 'request');

            ctrl.pollJobStatus();

            expect(Ext.Ajax.request).toHaveBeenCalledWith({
                url: '/status/3ad25048-26f6-11e1-851e-00012e260791.json',
                scope: ctrl,
                success: jasmine.any(Function),
                failure: jasmine.any(Function)
            });
        });

        it('callback_jobrunning', function() {
            spyOn(Ext.Ajax, 'request').andCallFake(function(options) {
                 Ext.callback(options.success, ctrl, [{responseText:"{\"status\":\"RUNNING\"}"}]);
            });
            spyOn(annotateActionButton, 'setTooltip');
            spyOn(ctrl, 'getAnnotateActionButton').andReturn(annotateActionButton);

            ctrl.pollJobStatus();

            expect(annotateActionButton.setTooltip).toHaveBeenCalledWith('Job RUNNING, waiting for completion');
        });

        it('callback_jobstopped', function() {
            spyOn(Ext.Ajax, 'request').andCallFake(function(options) {
                 Ext.callback(options.success, ctrl, [{responseText:"{\"status\":\"STOPPED\"}"}]);
            });
            ctrl.pollTask = 1234;
            spyOn(Ext.TaskManager, 'stop');
            spyOn(annotateActionButton, 'setTooltip');
            spyOn(annotateActionButton, 'setIconCls');
            spyOn(annotateActionButton, 'setText');
            spyOn(annotateActionButton, 'setHandler');
            spyOn(annotateActionButton, 'enable');
            spyOn(ctrl, 'getAnnotateActionButton').andReturn(annotateActionButton);

            ctrl.pollJobStatus();

            expect(Ext.TaskManager.stop).toHaveBeenCalledWith(1234);
            expect(ctrl.pollTask).toBeUndefined();
            expect(annotateActionButton.setTooltip).toHaveBeenCalledWith('Job completed, fetch results');
            expect(annotateActionButton.setIconCls).toHaveBeenCalledWith('');
            expect(annotateActionButton.setText).toHaveBeenCalledWith('Fetch result');
            expect(annotateActionButton.setHandler).toHaveBeenCalledWith(jasmine.any(Function));
            expect(annotateActionButton.enable).toHaveBeenCalledWith();
        });

        it('callback_failure', function() {
            spyOn(Ext.Ajax, 'request').andCallFake(function(options) {
                 Ext.callback(options.failure, ctrl);
            });
            ctrl.pollTask = 1234;
            spyOn(Ext.TaskManager, 'stop');
            spyOn(Ext.Error, 'raise');

            ctrl.pollJobStatus();

            expect(ctrl.pollTask).toBeUndefined();
            expect(Ext.TaskManager.stop).toHaveBeenCalledWith(1234);
            expect(Ext.Error.raise).toHaveBeenCalledWith('Failed to poll job status');
        });
     });

    describe('assign_struct2peakAction', function() {
        var button;

        beforeEach(function() {
            var params = {
                scanid: 1089,
                molid: 92,
                mz: 463.08856201171875
            };
            button = { pressed: true, params: params};
        });

        it('assign', function() {
            spyOn(Ext.Ajax, 'request');

            ctrl.assign_struct2peakAction(button);

            expect(Ext.Ajax.request).toHaveBeenCalledWith({
                url: '/rpc/3ad25048-26f6-11e1-851e-00012e260790/assign',
                params: button.params,
                success: jasmine.any(Function),
                failure: jasmine.any(Function)
            });
        });

        it('unassign', function() {
            spyOn(Ext.Ajax, 'request');
            button.pressed = false;

            ctrl.assign_struct2peakAction(button);

            expect(Ext.Ajax.request).toHaveBeenCalledWith({
                url: '/rpc/3ad25048-26f6-11e1-851e-00012e260790/unassign',
                params: button.params,
                success: jasmine.any(Function),
                failure: jasmine.any(Function)
            });
        });

        it('callback_success', function() {
            var f = { callback: function() {} };
            spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
            Ext.util.Observable.capture(ctrl.application, f.callback);

            spyOn(Ext.Ajax, 'request').andCallFake(function(options) {
                 Ext.callback(options.success);
            });

            ctrl.assign_struct2peakAction(button);

            expect(f.callback).toHaveBeenCalledWith('assignmentchanged', button.pressed, button.params);
            Ext.util.Observable.releaseCapture(ctrl.application);
        });

        it('callback_failure', function() {
            var f = { callback: function() {}};
            spyOn(f, 'callback').andReturn(false); // listeners dont hear any events
            Ext.util.Observable.capture(ctrl.application, f.callback);
            spyOn(Ext.Ajax, 'request').andCallFake(function(options) {
                 Ext.callback(options.failure);
            });
            spyOn(Ext.Error, 'raise');

            ctrl.assign_struct2peakAction(button);

            expect(Ext.Error.raise).toHaveBeenCalledWith('Failed to (un)assign molecule to peak');
            expect(f.callback).not.toHaveBeenCalledWith('assignmentchanged', button.pressed, button.params);

            Ext.util.Observable.releaseCapture(ctrl.application);
        });
     });

    it('onLaunch', function() {
    	spyOn(ctrl, 'applyRole');

    	ctrl.onLaunch();

    	expect(ctrl.applyRole).toHaveBeenCalledWith();
    });

    it('canassign', function() {
  	  ctrl.application.features.assign = true;
      assignbut = jasmine.createSpyObj('abut', [ 'hide']);
      spyOn(ctrl,'getAssignStruct2PeakButton').andReturn(assignbut);

  	  ctrl.applyRole();

  	  expect(assignbut.hide).not.toHaveBeenCalledWith();
    });

    it('cantassign', function() {
  	  ctrl.application.features.assign = false;
      assignbut = jasmine.createSpyObj('abut', [ 'hide']);
      spyOn(ctrl,'getAssignStruct2PeakButton').andReturn(assignbut);

  	  ctrl.applyRole();

  	  expect(assignbut.hide).toHaveBeenCalledWith();
    });


    it('showHelp', function() {
       spyOn(ctrl.application, 'showHelp');

       ctrl.showHelp();

       expect(ctrl.application.showHelp).toHaveBeenCalledWith('substructures');
    });

  });
});