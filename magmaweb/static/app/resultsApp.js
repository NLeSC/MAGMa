/**
 * MAGMaWeb results application
 * @author <a href="mailto:s.verhoeven@esciencecenter.nl">Stefan Verhoeven</a>
 *
 *
 * Example config:
 *
 *     @example
 *     app = Ext.create('Esc.magmaweb.resultsApp', {
 *       maxmslevel: 3,
 *       urls: {
 *         home: '/',
 *         fragments: '/fragments/{0}/{1}.json',
 *         mspectra: '/mspectra/{0}.json?mslevel={1}',
 *         extractedionchromatogram: '/extractedionchromatogram/{0}.json',
 *         metabolites: '/metabolites.json',
 *         chromatogram: '/chromatogram.json'
 *       }
 *     });
 *
 * Note! Example requires that Esc.magmaweb, Esc namespaces to be resolvable.
 */
Ext.define('Esc.magmaweb.resultsApp', {
  extend:'Ext.app.Application',
  constructor: function(config) {
    console.log('Construct app');
    this.initConfig(config);
    this.callParent(arguments);
    return this;
  },
  name: 'Esc.magmaweb',
  controllers: [ 'Metabolites', 'Fragments', 'Scans', 'MSpectras' ],
  requires: [
    'Ext.panel.Panel',
    'Ext.container.Viewport',
    'Ext.layout.container.Border',
    'Ext.Img',
    'Ext.window.Window',
    'Ext.form.Panel',
    'Esc.magmaweb.view.fragment.AnnotateFieldSet'
  ],
  config: {
    /**
     * Metabolite grid page size.
     * @cfg {Number}
     */
    pageSize: 10,
    /**
     * Maximum MS level or nr of MS levels.
     * @cfg {Number}
     */
    maxmslevel: 2,
    /**
     * Job identifier
     * @cfg String
     */
    jobid: null,
    /**
     * Endpoints/templates for contacting server.
     * @cfg {Object}
     */
    urls: {
        /**
         * Homepage.
         * @cfg {String} urls.home
         */
        home: null,
        /**
         * Fragments endpoint.
         * Tokenized string with scanid and metid tokens.
         * @cfg {String} urls.fragments
         */
        fragments: null,
        /**
         * MSpectra endpoint.
         * Tokenized string with scanid and mslevel tokens.
         * @cfg {String} urls.mspectra
         */
        mspectra: null,
        /**
         * Extracted ion chromatogram endpoint.
         * Tokenized string with metid token.
         * @cfg {String} urls.extractedionchromatogram
         */
        extractedionchromatogram: null,
        /**
         * Metabolites endpoint.
         * @cfg {String} urls.metabolites
         */
        metabolites: null,
        /**
         * Chromatogram endpoint.
         * @cfg {String} urls.chromatogram
         */
        chromatogram: null,
        /**
         * Stderr endpoint.
         * @cfg {String} urls.stderr
         */
        stderr: null
    }
  },
  /**
   * when a metabolite and scan are selected then load fragments
   * @property {Object} selected
   * @property {Boolean} selected.scanid Scan identifier
   * @property {Boolean} selected.metid Metabolite identifier
   */
  selected: { scanid: false, metid: false },
  /**
   * Can only annotate when there are structures and ms data.
   * @property {Object} annotatable
   * @property {Boolean} annotabable.structures Whether there are structures
   * @property {Boolean} annotabable.msdata Whether there is ms data
   */
  annotatable: { structures: false, msdata: false },
  /**
   * Shows annotate form in modal window
   */
  showAnnotateForm: function() {
    var me = this;
    if (!this.annotateForm) {
        this.annotateForm = Ext.create('Ext.window.Window', {
            title: 'Annotate all structures',
            modal: true,
            height: 400,
            width: 600,
            layout: 'fit',
            closeAction: 'hide',
            items: {
                xtype: 'form',
                bodyPadding: 5,
                defaults: { bodyPadding: 5 },
                border: false,
                autoScroll: true,
                url: me.rpcUrl('annotate'),
                items: [{
                    xtype : 'annotatefieldset'
                }],
                buttons: [{
                  text: 'Submit',
                  handler: function() {
                    var form = this.up('form').getForm();
                    var wf = this.up('window');
                    if (form.isValid()) {
                      form.submit({
                        waitMsg: 'Submitting action ...',
                        success: function(fp, o) {
                          var response = Ext.JSON.decode(o.response.responseText);
                          me.fireEvent('rpcsubmitsuccess', response.jobid);
                          wf.hide();
                        },
                        failure: function(form, action) {
                          console.log(action.failureType);
                          console.log(action.result);
                          wf.hide();
                        }
                      });
                    }
                  }
                }, {
                  text: 'Reset',
                  handler: function() {
                    this.up('form').getForm().reset();
                  }
                }]
            }
        });
    }
    this.annotateForm.show();
  },
  /**
   * Logs error in console and shows a error message box to user
   *
   * @param {Ext.Error} err The raised error
   */
  errorHandle: function(err) {
      console.error(err);
      Ext.Msg.show({
          title: 'Error',
          msg: err.msg,
          buttons: Ext.Msg.OK,
          icon: Ext.Msg.ERROR
      });
      return true;
  },
  /**
   * Get url of rpc method
   * @param {String} method
   * @return {Url}
   */
  rpcUrl: function(method) {
    return this.urls.home+'rpc/'+this.jobid+'/'+method;
  },
  /**
   * Pols status of jobid and when completed redirects to results of new job.
   * @param {String} jobid Identifier of new job
   */
  rpcSubmitted: function(jobid) {
      var me = this;
      me.newjobid = jobid;
      // Overwrite annotate button to waiting/cancel button
      var annot_button = Ext.getCmp('annotateaction');
      annot_button.setIconCls('icon-loading');
      annot_button.setText('Waiting');
      annot_button.setTooltip('Job submitted, waiting for completion');
      annot_button.setHandler(function() {
        Ext.MessageBox.confirm('Cancel job', 'Job is still running. Do you want to cancel it?', function(but) {
            if (but == 'yes') {
                // TODO cancel job
                console.log('Cancelling job');
            }
        }, me);
      });
      // Start polling
      me.pollTask = Ext.TaskManager.start({
        run: function() {
            Ext.Ajax.request({
              url: me.urls.home+'status/'+me.newjobid+'.json',
              success: function(o) {
                var response = Ext.JSON.decode(o.responseText);
                console.log(response.status);
                if (response.status == 'STOPPED') {
                  Ext.TaskManager.stop(me.pollTask);
                  delete me.pollTask;
                  annot_button.setIconCls('');
                  annot_button.setText('Fetch result');
                  annot_button.setTooltip('Job completed, fetch results');
                  annot_button.setHandler(function() {
                    Ext.MessageBox.confirm('Fetch result', 'Job has been completed. Do you want to fetch results?', function(but) {
                        if (but == 'yes') {
                            window.location = me.urls.home+'results/'+me.newjobid;
                        }
                    }, me);
                  });
                } else {
                    annot_button.setTooltip('Job '+response.status+', waiting for completion');
                }
              },
              failure: function() {
                Ext.TaskManager.stop(me.pollTask);
                delete me.pollTask;
              },
              scope: me
            });
        },
        interval: 5000
      });
  },
  /**
   * Creates mspectraspanels and viewport and fires/listens for mspectra events
   * Registers error handle
   */
  launch: function() {
    Ext.Error.handle = this.errorHandle;
    var me = this;
    this.addEvents(
      /**
       * @event
       * Triggered when a metabolite and scan are selected together.
       * @param {Number} scanid Scan identifier.
       * @param {Number} metid Metabolite identifier.
       */
      'scanandmetaboliteselect',
      /**
       * @event
       * Triggered when a metabolite and scan are no longer selected together.
       * @param {Number} scanid Scan identifier.
       * @param {Number} metid Metabolite identifier.
       */
      'scanandmetabolitenoselect',
      /**
       * @event
       * Triggered when a rpc method has been submitted successfully.
       * @param {String} jobid Job identifier of new job submitted.
       */
      'rpcsubmitsuccess'
    );

    // uncomment to see all application events fired in console
//    Ext.util.Observable.capture(this, function() { console.log(arguments);return true;});

    this.on('metaboliteselect', function(metid) {
      this.selected.metid = metid;
      if (this.selected.metid && this.selected.scanid) {
        this.fireEvent('scanandmetaboliteselect', this.selected.scanid, metid);
      }
    }, this);
    this.on('selectscan', function(scanid) {
        this.selected.scanid = scanid;
        if (this.selected.metid && this.selected.scanid) {
            this.fireEvent('scanandmetaboliteselect', scanid, this.selected.metid);
        }
    }, this);
    this.on('noselectscan', function() {
        if (this.selected.metid && this.selected.scanid) {
            this.fireEvent('scanandmetabolitenoselect');
        }
        this.selected.scanid = false;
    }, this);
    this.on('metabolitedeselect', function() {
        if (this.selected.metid && this.selected.scanid) {
            this.fireEvent('scanandmetabolitenoselect');
        }
        this.selected.metid = false;
    }, this);
    this.on('metabolitenoselect', function() {
        if (this.selected.metid && this.selected.scanid) {
            this.fireEvent('scanandmetabolitenoselect');
        }
        this.selected.metid = false;
    }, this);

    this.on('mspectraload', function(scanid, mslevel) {
      Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Level '+mslevel+' scan '+scanid);
    });
    this.on('mspectraclear', function(mslevel) {
      Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Level '+mslevel+' scan ...');
    });
    this.on('peakmouseover', function(peak, mslevel, scanid) {
      Ext.getCmp('mspectra'+mslevel+'panel').header.setTitle('Level '+mslevel+' scan '+scanid+' (m/z='+peak.mz+', intensity='+peak.intensity+')');
    });

    this.on('metaboliteload', function(store) {
        this.annotatable.structures = store.getTotalUnfilteredCount() > 0;
        if (this.annotatable.structures && this.annotatable.msdata) {
            Ext.getCmp('annotateaction').enable();
        } else {
            Ext.getCmp('annotateaction').disable();
        }
    }, this);
    this.on('chromatogramload', function(chromatogram) {
        this.annotatable.msdata = chromatogram.data.length > 0;
        if (this.annotatable.structures && this.annotatable.msdata) {
            Ext.getCmp('annotateaction').enable();
        } else {
            Ext.getCmp('annotateaction').disable();
        }
    }, this);
    this.on('rpcsubmitsuccess', this.rpcSubmitted, this);

    console.log('Launch app');

    var msspectrapanels = [];
    var mspectras = this.getController('MSpectras').mspectras;
    for (var mslevel = 1; mslevel <= this.getMaxmslevel(); mslevel++) {
      msspectrapanels.push({
        title: 'Level '+mslevel+' scan ...',
        id: 'mspectra'+mslevel+'panel',
        collapsible: true,
        tools: [{
          type: 'restore',
          tooltip: 'Center level '+mslevel+' scan',
          disabled: true,
          action: 'center'
        }],
        items: mspectras[mslevel]
      });
    }
    if (this.getMaxmslevel() > 0) {
        var mspectrapanel = Ext.create('Ext.panel.Panel', {
            region:'south',
            height: '50%',
            split: true,
            collapsible: true,
            hideCollapseTool: true,
            border: false,
            preventHeader: true,
            id: 'mspectrapanel',
            layout: {
              type: 'vbox',
              align: 'stretch'
            },
            defaults: {
              flex: 1,
              layout:'fit',
              border: false
            },
            items: msspectrapanels
        });
    } else {
        var mspectrapanel = Ext.create('Ext.panel.Panel', {
            region:'south',
            height: '50%',
            split: true,
            collapsible: true,
            hideCollapseTool: true,
            border: false,
            title: 'Scans',
            id: 'mspectrapanel',
            html: 'No scans available: Upload ms data'
        });
    }

    var infoWindow = Ext.create('Ext.window.Window', {
        title: 'Information',
        width: 600,
        autoHeight: true,
        closeAction: 'hide',
        contentEl: 'resultsinfo',
        tools: [{
            type: 'gear',
            tooltip: 'Set description',
            handler: function() {
                Ext.MessageBox.prompt('Description', 'Please enter a description:', function(button, description) {
                    if (button == 'ok') {
                        Ext.Ajax.request({
                            url: me.rpcUrl('set_description'),
                            params: {
                                description: description
                            },
                            success: function(response) {
                                Ext.get('description').setHTML(description);
                            }
                        });
                    }
                }, this, true, Ext.get('description').getHTML());
            }
        }]
    });

    // header
    var header_side = Ext.create('Ext.panel.Panel', {
      region: 'north',
      layout: {
        type: 'hbox',
        align: 'middle',
        padding: 2
      },
      items: [{
        xtype: 'component',
        cls: 'x-logo',
        html: '<a href="'+me.urls.home+'" data-qtip="<b>M</b>s <b>A</b>nnotation based on in silico <b>G</b>enerated <b>M</b>et<b>a</b>bolites">MAGMa</a>'
      }, {
        xtype:'tbspacer',
        flex:1 // aligns buttongroup right
      }, {
          xtype: 'buttongroup',
          columns: 3,
          items: [{
              text: 'Restart',
              handler: function() {
                window.location = me.urls.home;
              },
              tooltip: 'Upload a new dataset'
            },{
              text: 'Download',
              tooltip: 'Download the different results files',
              menu: {
                  items: [{
                      text: 'Metabolites',
                      handler: function() {
                          // TODO replace handler with action
                          me.getController('Metabolites').download();
                      }
                  }, {
                      text: 'Fragments',
                      disabled: true
                  }, {
                      text: 'Error log',
                      href: me.urls.stderr,
                      hrefTarget: '_new'
                  }]
              }
            },{
              text: 'Annotate',
              tooltip: 'Annotate all structures',
              id: 'annotateaction',
              iconCls: 'icon-annot',
              handler: me.showAnnotateForm.bind(me),
              disabled: true
            },{
              text: 'Help',
              tooltip: 'Goto help pages',
              disabled: true
            }, {
              text: 'Information',
              tooltip: 'Information about analysis parameters',
              handler: function() {
                  infoWindow.show();
              }
            }]
        }]
    });

    var master_side = Ext.create('Ext.panel.Panel', {
      // master side
      region: 'center',
      layout: 'border',
      border: false,
      items:[{
        region:'center',
        border: false,
        xtype: 'metabolitelist'
      },{
        region:'south',
        hideCollapseTool: true,
        collapsible: true,
        height: '50%',
        split: true,
        xtype: 'chromatogrampanel',
        border: false
      }]
    });

    // detail side
    var detail_side = Ext.create('Ext.panel.Panel', {
      region: 'east',
      split: true,
      collapsible: true,
      layout: 'border',
      width: 600,
      hideCollapseTool: true,
      border: false,
      preventHeader: true,
      items:[{
        region: 'center',
        xtype: 'fragmenttree',
        border: false
      },
      mspectrapanel
      ]
    });

    Ext.create('Ext.Viewport', {
      layout: 'border',
      items:[ master_side, detail_side, header_side ]
    });
  }
});
