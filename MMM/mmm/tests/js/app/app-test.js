Ext.Loader.setConfig({
  enabled: true,
//  disableCaching: false, // otherwise can not use firebug breakpoint
  paths: {
    'Esc.mmm': '../../../static/app',
    'Esc': '../../../static/esc',
    'Ext.ux': '../../../static/ext-4.0.7/examples/ux'
  }
});

var Application = null;

Ext.onReady(function() {
  Application = Ext.create('Esc.mmm.resultsApp', {
      appFolder: "../../../static/app",
      maxmslevel: 3,
      ms_intensity_cutoff: 2000000,
      urls: {
        home: '/',
        nlesclogo: '/',
        fragments: 'data/fragments.s{0}.m{1}.json',
        mspectra: 'data/mspectra.{0}.json?mslevel={1}',
        extractedionchromatogram: 'data/extractedionchromatogram.{0}.json',
        metabolites: 'data/metabolites.json',
        chromatogram: 'data/chromatogram.json'
      },
      launch: function() {
          //include the tests in the test.html head
          jasmine.getEnv().addReporter(new jasmine.TrivialReporter());
          jasmine.getEnv().execute();
      }
  });
});
