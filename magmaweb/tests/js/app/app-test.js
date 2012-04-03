Ext.Loader.setConfig({
  enabled: true,
//  disableCaching: false, // otherwise can not use firebug breakpoint
  paths: {
    'Esc.magmaweb': '../../../static/app',
    'Esc': '../../../static/esc',
    'Ext.ux': '../../../static/extjs-4.1.0-rc1/examples/ux'
  }
});

var Application = null;

Ext.onReady(function() {
  Application = Ext.create('Esc.magmaweb.resultsApp', {
      appFolder: "../../../static/app",
      maxmslevel: 3,
      ms_intensity_cutoff: 2000000,
      urls: {
        home: '/',
        fragments: 'data/fragments.s{0}.m{1}.json',
        mspectra: 'data/mspectra.{0}.json?mslevel={1}',
        extractedionchromatogram: 'data/extractedionchromatogram.{0}.json',
        metabolites: 'data/metabolites.json',
        metabolitescsv: 'data/metabolites.csv',
        chromatogram: 'data/chromatogram.json',
        stderr: 'data/stderr.txt'
      },
      launch: function() {
          //include the tests in the test.html head
          jasmine.getEnv().addReporter(new jasmine.TrivialReporter());
          jasmine.getEnv().execute();
      }
  });
});
