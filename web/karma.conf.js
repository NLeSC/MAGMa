// Karma configuration
// Generated on Fri Aug 23 2013 17:04:50 GMT+0200 (CEST)

module.exports = function(config) {
  config.set({

    // base path, that will be used to resolve files and exclude
    basePath: '',


    // frameworks to use
    frameworks: ['jasmine'],


    // list of files / patterns to load in the browser
    files: [
      // all js files used in index.html
      'magmaweb/static/extjs-4.2.0/ext-all-debug.js',
      'magmaweb/static/extjs-4.2.0/examples/ux/grid/filter/Filter.js',
      'magmaweb/static/extjs-4.2.0/examples/ux/grid/filter/DateFilter.js',
      'magmaweb/static/extjs-4.2.0/examples/ux/grid/**/*.js',
      'magmaweb/static/ChemDoodleWeb/ChemDoodleWeb-libs.js',
      'magmaweb/static/ChemDoodleWeb/ChemDoodleWeb.js',
      'magmaweb/static/d3/d3.min.js',
      // app
      'magmaweb/tests/js/app/loader.js',
      'magmaweb/static/esc/**/*.js',
      'magmaweb/static/app/model/*.js',
      'magmaweb/static/app/store/*.js',
      'magmaweb/static/app/view/scan/UploadFieldSet.js',
      'magmaweb/static/app/view/metabolite/ReactionFilter.js',
      'magmaweb/static/app/view/metabolite/ReactionColumn.js',
      'magmaweb/static/app/view/metabolite/ScenarioField.js',
      'magmaweb/static/app/view/metabolite/MetabolizeFieldSet.js',
      'magmaweb/static/app/view/fragment/AnnotateFieldSet.js',
      'magmaweb/static/app/view/scan/UploadForm.js',
      'magmaweb/static/app/view/*/*.js',
      'magmaweb/static/app/view/Viewport.js',
      'magmaweb/static/app/controller/*.js',
      'magmaweb/static/app/resultsApp.js',
      'magmaweb/tests/js/app/app-test.js',
      // tests
      {pattern: 'magmaweb/tests/js/app/data/*.*', watched: false, included: false, served: true},
      'magmaweb/tests/js/specs/*.spec.js',
      'magmaweb/tests/js/app/specs/*.js'
    ],

    preprocessors: {
        // source files, that you wanna generate coverage for
        // do not include tests or libraries
        // (these files will be instrumented by Istanbul)
//        'magmaweb/static/esc/**/*.js': ['coverage'],
//        'magmaweb/static/app/**/*.js': ['coverage'],
//        'magmaweb/static/app/resultsApp.js': ['coverage']
    },

    // list of files to exclude
    exclude: [

    ],

    // test results reporter to use
    // possible values: 'dots', 'progress', 'junit', 'growl', 'coverage'
    reporters: ['dots', 'junit', 'coverage'],

    junitReporter: {
        outputFile: 'reports/TEST-results.xml'
    },
    coverageReporter: {
        dir: 'coverage/',
        reporters: [{
            type: 'lcov' // for viewing html pages and SonarQube
        }, {
            type: 'cobertura' // for use in Jenkins
        }]
    },

    // web server port
    port: 9876,


    // enable / disable colors in the output (reporters and logs)
    colors: true,


    // level of logging
    // possible values: config.LOG_DISABLE || config.LOG_ERROR || config.LOG_WARN || config.LOG_INFO || config.LOG_DEBUG
    logLevel: config.LOG_INFO,


    // enable / disable watching file and executing tests whenever any file changes
    autoWatch: false,


    // Start these browsers, currently available:
    // - Chrome
    // - ChromeCanary
    // - Firefox
    // - Opera
    // - Safari (only Mac)
    // - PhantomJS
    // - IE (only Windows)
    browsers: ['PhantomJS'],


    // If browser does not capture in given timeout [ms], kill it
    captureTimeout: 60000,


    // Continuous Integration mode
    // if true, it capture browsers, run tests and exit
    singleRun: true
  });
};
