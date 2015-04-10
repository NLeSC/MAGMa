Ext.require('Esc.magmaweb.resultsApp');

var Application = null;

Ext.onReady(function() {
  var appRootBase = 'base/magmaweb/tests/js/app';

  Application = Ext.create('Esc.magmaweb.resultsApp', {
    maxmslevel: 3,
    jobid: '3ad25048-26f6-11e1-851e-00012e260790',
    urls: {
      home: appRootBase,
      fragments: appRootBase + '/data/fragments.s{0}.m{1}.json',
      mspectra: appRootBase + '/data/mspectra.{0}.json?mslevel={1}',
      extractedionchromatogram: appRootBase + '/data/extractedionchromatogram.{0}.json',
      chromatogram: appRootBase + '/data/chromatogram.json',
      stderr: appRootBase + '/data/stderr.txt'
    },
    launch: function() {
      // do nothing, karma will start the tests for us
    },
    // mock runinfo url to a static file
    runInfoUrl: function() {
      return appRootBase + '/data/runinfo.json';
    },
    // use test files for molecules in data/
    moleculesUrl: function(format) {
      return appRootBase + 'data/molecules.' + format;
    },
    showHelp: function(section) {}
  });

  Ext.getBody().createChild({
    tag: 'div',
    id: 'resultsinfo',
    html: 'placeholder for resultsinfo'
  });
  Ext.getBody().createChild({
    tag: 'div',
    id: 'logos',
    html: 'placeholder for logos'
  });
});
