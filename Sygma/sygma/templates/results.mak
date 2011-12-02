<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <title>MSygma - Results</title>
  <meta http-equiv="Content-Type" content="text/html;charset=UTF-8"/>
<link rel="stylesheet" href="${request.static_url('sygma:static/ChemDoodleWeb/ChemDoodleWeb.css')}" type="text/css"></link>
<link rel="stylesheet" href="${request.static_url('sygma:static/ext-4.0.7-gpl/resources/css/ext-all.css')}" type="text/css"></link>
<link rel="stylesheet" type="text/css" href="${request.static_url('sygma:static/ext-4.0.7-gpl/examples/ux/grid/css/GridFilters.css')}" />
<link rel="stylesheet" type="text/css" href="${request.static_url('sygma:static/ext-4.0.7-gpl/examples/ux/grid/css/RangeMenu.css')}" />
<script type="text/javascript" src="${request.static_url('sygma:static/ext-4.0.7-gpl/bootstrap.js')}"></script>
<script type="text/javascript" src="${request.static_url('sygma:static/ChemDoodleWeb/ChemDoodleWeb-libs.js')}"></script>
<script type="text/javascript" src="${request.static_url('sygma:static/ChemDoodleWeb/ChemDoodleWeb.js')}"></script>
<script type="text/javascript" src="${request.static_url('sygma:static/d3/d3.js')}"></script>
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

</style><script type="text/javascript">

Ext.Loader.setConfig({
  enabled: true,
//  disableCaching: false, // uncomment to use firebug breakpoints
  paths: {
    'Esc.msygma': '${request.static_url('sygma:static/app')}',
    'Esc': '${request.static_url('sygma:static/esc')}',
    'Ext.ux': '${request.static_url('sygma:static/ext-4.0.7-gpl/examples/ux')}'
  }
});

app = Ext.create('Esc.msygma.resultsApp', {
  appFolder: "${request.static_url('sygma:static/app')}",
  maxmslevel: ${maxmslevel},
  ms_intensity_cutoff: ${run.ms_intensity_cutoff},
  urls: {
    fragments: '${request.application_url}/fragments/{0}/{1}.json',
    mspectra: '${request.application_url}/mspectra/{0}.json?mslevel={1}',
    extractedionchromatogram: '${request.application_url}/extractedionchromatogram/{0}.json',
    metabolites: '${request.route_url('metabolites.json')}',
    chromatogram: '${request.route_url('chromatogram.json')}'
  }
});
</script>
</head>
<body>
</body>
</html>

