Ext.Loader.setConfig({
  enabled: true,
//  disableCaching: false, // otherwise can not use firebug breakpoint
  paths: {
    'Esc.magmaweb': '../../../static/app',
    'Esc': '../../../static/esc',
    'Ext.ux': '../../../static/extjs-4.2.0/examples/ux'
  }
});

Ext.require('Esc.magmaweb.resultsApp');

var Application = null;

Ext.onReady(function() {
  Application = Ext.create('Esc.magmaweb.resultsApp', {
      maxmslevel: 3,
      jobid: '3ad25048-26f6-11e1-851e-00012e260790',
      urls: {
        home: '/',
        fragments: 'data/fragments.s{0}.m{1}.json',
        mspectra: 'data/mspectra.{0}.json?mslevel={1}',
        extractedionchromatogram: 'data/extractedionchromatogram.{0}.json',
        chromatogram: 'data/chromatogram.json',
        stderr: 'data/stderr.txt'
      },
      launch: function() {
          //include the tests in the test.html head
          jasmine.getEnv().addReporter(new jasmine.TrivialReporter());
          jasmine.getEnv().addReporter(new jasmine.JUnitXmlReporter());
          jasmine.getEnv().execute();
      },
      // mock runinfo url to a static file
      runInfoUrl: function() { return 'data/runinfo.json' },
      // use test files for molecules in data/
      moleculesUrl: function(format) {
          return 'data/molecules.'+format;
      },
      showHelp: function(section) {}
  });
});
