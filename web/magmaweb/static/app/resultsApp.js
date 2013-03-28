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
  name: 'Esc.magmaweb',
  appFolder: Ext.Loader.getPath('Esc.magmaweb'),
  extend:'Ext.app.Application',
  constructor: function(config) {
    Ext.log({}, 'Construct app');
    this.initConfig(config);
    this.callParent(arguments);
    return this;
  },
  autoCreateViewport: true,
  controllers: [ 'Metabolites', 'Fragments', 'Scans', 'MSpectras' ],
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
     * Feature toggles.
     * @cfg {Object}
     */
    features: {
        /**
         * Whether user has rights to (un)assign a peak to a structure.
         * @cfg {Boolean}
         */
        run: false,
        /**
         * Whether user has rights to (un)assign a peak to a structure.
         * @cfg {Boolean}
         */
        assign: true,
        /**
         * Whether user is logged in, shows login or logout button.
         * @cfg {Boolean}
         */
        authenticated: true,
        /**
         * Whether user is needs to be logged in, shows login or logout button.
         * @cfg {Boolean}
         */
        anonymous: false,
        /**
         * Whether restrictions should be applied ie force one spectral tree
         * @cfg {Boolean}
         */
        restricted: false
    },
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
      Ext.log({level:'error'}, err);
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
   * @return {String}
   */
  rpcUrl: function(method) {
    return this.urls.home+'rpc/'+this.jobid+'/'+method;
  },
  /**
   * Get url of runinfo json, used to set defaults in forms.
   * @return {String}
   */
  runInfoUrl: function() {
    return this.urls.home+'results/'+this.jobid+'/runinfo.json';
  },
  /**
   * Get metabolites url based on format.
   *
   * @param {String} format Can be json, csv or sdf.
   * @return {String}
   */
  metabolitesUrl: function(format) {
    return this.urls.home+'results/'+this.jobid+'/metabolites.'+format;
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

    Ext.log({}, 'Launch app');

    /**
     * @property {Ext.window.Window} infoWindow
     * Information window which shows settings used, description and error log.
     */
    this.infoWindow = Ext.create('Ext.window.Window', {
        title: 'Information',
        width: 600,
        autoHeight: true,
        closeAction: 'hide',
        contentEl: 'resultsinfo',
        tools: [{
            type: 'save',
            tooltip: 'Save log file',
            handler: function() {
                window.open(me.urls.stderr, 'Log');
            }
        }]
    });
    // cant use this.control, find component and setHandler
    Ext.ComponentQuery.query('component[action=information]')[0].setHandler(function() {
        me.infoWindow.show();
    });
  },
  /**
   * Apply role to user interface.
   * Checks canRun and if false removes all action buttons.
   * All controllers should apply roles to user themselves if required.
   */
  applyRole: function() {
    var features = this.features;
    if (!this.features.run) {
        Ext.ComponentQuery.query('component[id=annotateaction]')[0].hide();
    }
    if (features.authenticated || features.anonymous) {
        Ext.ComponentQuery.query('component[text=Login]')[0].hide();
        if (features.anonymous) {
          // non-anonymous authenticated can not logout
          Ext.ComponentQuery.query('component[text=Logout]')[0].hide();
        }
    } else {
        Ext.ComponentQuery.query('component[text=Logout]')[0].hide();
        Ext.ComponentQuery.query('component[text=Workspace]')[0].hide();
    }
  }
});
