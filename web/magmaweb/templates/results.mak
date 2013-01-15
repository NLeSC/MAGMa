<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Results</title>
<link rel="stylesheet"
	href="${request.extjsroot}/resources/css/ext-all.css" type="text/css"></link>
<link rel="stylesheet"
	href="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.9.2.custom.css')}"
	type="text/css"></link>
<link rel="stylesheet" type="text/css"
	href="${request.extjsroot}/examples/ux/grid/css/GridFilters.css" />
<link rel="stylesheet" type="text/css"
	href="${request.extjsroot}/examples/ux/grid/css/RangeMenu.css" />
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

.y.axis line,.y.axis path,.x.axis line,.x.axis path {
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

.marker {
	stroke: lightgray;
	fill: none;
	stroke-width: 1.5px;
	cursor: pointer;
}

.annotated {
	stroke: lightgreen;
}

line.assigned {
	stroke: orange !important; /* #afdddd; */
}

.annotatedandassigned {
	stroke: #F4A460;
}

.selected {
	fill: darkgreen !important;
}

line.mspeak {
	stroke-width: 1px;
	stroke: black;
}

.fragmenttree .x-grid-cell-inner {
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
	text-decoration: none;
	margin-left: auto;
	margin-right: auto;
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
	background-image:
		url(${request.extjsroot}/examples/shared/icons/fam/cog.png);
}

.icon-add {
	background-image:
		url(${request.extjsroot}/examples/writer/images/add.png) !important;
}

.icon-annot {
	background-image:
		url(${request.extjsroot}/examples/feed-viewer/images/post_go.gif);
}

.icon-loading {
	background-image:
		url(${request.extjsroot}/resources/themes/images/default/grid/loading.gif);
	cursor: wait;
}

.icon-connect {
	background-image:
		url(${request.extjsroot}/examples/shared/icons/fam/connect.gif)
		!important;
}
</style>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb-libs.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/ChemDoodleWeb.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/jquery-ui-1.9.2.custom.min.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/ChemDoodleWeb/sketcher/ChemDoodleWeb-sketcher.js')}"></script>
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/d3/d3.v3.min.js')}"></script>
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
## Comment out below for development or when running sencha build, every
## class is loaded when commented out
<script type="text/javascript"
	src="${request.static_url('magmaweb:static/app/resultsApp-all.js')}"></script>
<script type="text/javascript">

Ext.require('Esc.magmaweb.resultsApp');
<%!
import json
%>
var app;
Ext.onReady(function() {
    app = Ext.create('Esc.magmaweb.resultsApp', {
      maxmslevel: ${maxmslevel},
      jobid: '${jobid}',
      canRun: ${json.dumps(canRun)|n},
      is_user_authenticated: ${json.dumps(request.user is not None)},
      urls: {
        home: '${request.route_path('home')}',
        fragments: '${request.application_url}/results/${jobid}/fragments/{0}/{1}.json',
        mspectra: '${request.application_url}/results/${jobid}/mspectra/{0}.json?mslevel={1}',
        extractedionchromatogram: '${request.application_url}/results/${jobid}/extractedionchromatogram/{0}.json',
        chromatogram: '${request.route_path('chromatogram.json',jobid=jobid)}',
        stderr: '${request.route_path('stderr.txt',jobid=jobid)}'
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
		        'sketcher_canvas', 500, 300, {
		        	useServices: false, oneMolecule: true
		        }
			);
			sketcher.repaint();
			sketcher.toolbarManager.setup();
			sketcher.toolbarManager.buttonSave.disable();
			sketcher.toolbarManager.buttonOpen.disable();
		</script>
	</div>
</body>
</html>
