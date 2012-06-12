<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <title>MAGMa - Results</title>
  <meta http-equiv="Content-Type" content="text/html;charset=UTF-8"/>
<link rel="stylesheet" href="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb.css')}" type="text/css"></link>
<link rel="stylesheet" href="${request.extjsroot}/resources/css/ext-all.css" type="text/css"></link>
<link rel="stylesheet" href="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.8.7.custom.css')}" type="text/css"></link>
<link rel="stylesheet" type="text/css" href="${request.extjsroot}/examples/ux/grid/css/GridFilters.css" />
<link rel="stylesheet" type="text/css" href="${request.extjsroot}/examples/ux/grid/css/RangeMenu.css" />
<style type="text/css">

path.line {
  fill: none;
  stroke: #66d;
  stroke-width: 2px;
}

path.metaboliteline {
  fill: none;
  stroke: #6d6;
  stroke-width: 2px;
}

.axis {
  shape-rendering: crispEdges;
}

.y.axis line, .y.axis path, .x.axis line, .x.axis path {
  fill: none;
  stroke: #ccc;
}

.peaks {
  fill: none;
  stroke-width: 1px;
  stroke: "#111";
}

.cutoffline {
  fill: none;
  stroke-width: 1px;
  stroke: #ddd;
  stroke-dasharray: 5px, 5px;
}

line.peak {
  stroke-width: 2px;
  stroke: #eee;
}

path.marker {
  stroke: lightgreen;
  fill: none;
  stroke-width: 1.5px;
  cursor: pointer;
}

.selectedscan, .selectedpeak {
  stroke: darkgreen !important;
  fill: darkgreen !important;
}

line.mspeak {
  stroke-width: 1px;
  stroke: black;
}

.fragmenttree .x-grid-cell-inner{
  height: 106px !important;
}

.x-logo a {
  font-size: 40px;
  padding-left: 520px;
  padding-top: 3px; /* aligns app title with text in logo  */
  background: url(${request.static_url('magmaweb:static/ESCIENCE_log_B_nl_long_cyanblack.jpg')}) no-repeat 5px 4px;
}

.x-title a {
  font-size: 40px;
  font-weight: bold;
  color: #333;
  text-decoration:none;
  padding-left: 520px;
  padding-top: 3px; /* aligns app title with text in logo  */
}

#resultsinfo {
  padding: 5px;
}

.infotable td {
  padding-right: 30px;
}

svg {
  background: white;
}

/* shared styles for the ActionColumn icons */
.x-action-col-cell img {
    height: 16px;
    width: 16px;
    cursor: pointer;
}
.x-action-col-cell img.metabolize-col {
    background-image: url(${request.extjsroot}/examples/shared/icons/fam/cog.png);
}

.icon-add {
    background-image: url(${request.extjsroot}/examples/writer/images/add.png) !important;
}

.icon-annot {
    background-image: url(${request.extjsroot}/examples/feed-viewer/images/post_go.gif);
}

.icon-loading {
    background-image: url(${request.extjsroot}/resources/themes/images/default/grid/loading.gif);
    cursor: wait;
}

</style>
<script type="text/javascript" src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb-libs.js')}"></script>
<script type="text/javascript" src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb.js')}"></script>
<script type="text/javascript" src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.8.7.custom.min.js')}"></script>
<script type="text/javascript" src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/ChemDoodleWeb-sketcher.js')}"></script>
<script type="text/javascript" src="${request.static_url('magmaweb:static/d3/d3.v2.min.js')}"></script>
<script type="text/javascript" src="${request.extjsroot}/ext.js"></script>
<script type="text/javascript">
if (!window.console) window.console = {};
if (!window.console.log) window.console.log = function() {};

Ext.Loader.setConfig({
  enabled: true,
  //disableCaching: false, // uncomment to use firebug breakpoints
  paths: {
    'Esc.magmaweb': '${request.static_url('magmaweb:static/app')}',
    'Esc': '${request.static_url('magmaweb:static/esc')}',
    'Ext.ux': '${request.extjsroot}/examples/ux'
  }
});

</script>
## Comment out below for development or when running sencha build, every class is loaded when commented out
<script type="text/javascript" src="${request.static_url('magmaweb:static/app/resultsApp-all.js')}"></script>
<script type="text/javascript">

Ext.require('Esc.magmaweb.resultsApp');

Ext.onReady(function() {
    Ext.create('Esc.magmaweb.resultsApp', {
      appFolder: "${request.static_url('magmaweb:static/app')}",
      maxmslevel: ${maxmslevel},
      jobid: '${jobid}',
      urls: {
        home: '${request.route_url('home')}',
        fragments: '${request.application_url}/results/${jobid}/fragments/{0}/{1}.json',
        mspectra: '${request.application_url}/results/${jobid}/mspectra/{0}.json?mslevel={1}',
        extractedionchromatogram: '${request.application_url}/results/${jobid}/extractedionchromatogram/{0}.json',
        metabolites: '${request.route_url('metabolites.json',jobid=jobid)}',
        metabolitescsv: '${request.route_url('metabolites.csv',jobid=jobid)}',
        chromatogram: '${request.route_url('chromatogram.json',jobid=jobid)}',
        stderr: '${request.route_url('stderr.txt',jobid=jobid)}'
      }
    });
});

</script>
</head>
<body>
<%include file="runinfo.mak"/>
<div id="sketcher_content" class="x-hidden">
<script language="javascript">
var sketcher = new ChemDoodle.SketcherCanvas(
        'sketcher_canvas', 500, 300,
        '${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/icons/')}',
        ChemDoodle.featureDetection.supports_touch(), false);
sketcher.repaint();
sketcher.toolbarManager.buttonSave.disable();
sketcher.toolbarManager.buttonOpen.disable();
</script>
</div>
</body>
</html>
