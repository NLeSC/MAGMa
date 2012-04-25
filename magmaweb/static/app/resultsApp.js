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
    'Ext.toolbar.Spacer',
    'Ext.container.ButtonGroup'
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
   * Get url of runinfo json, used to set defaults in forms.
   * @return {Url}
   */
  runInfoUrl: function() {
    return this.urls.home+'results/'+this.jobid+'/runinfo.json'
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
        xtype: 'metabolitepanel'
      },{
        region:'south',
        hideCollapseTool: true,
        collapsible: true,
        height: '50%',
        split: true,
        xtype: 'scanpanel',
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
