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
<link rel="stylesheet"  href="${request.extjsroot}/examples/ux/css/ItemSelector.css" type="text/css"></link>
<link rel="stylesheet"  href="${request.extjsroot}/examples/writer/writer.css" type="text/css"></link>
<link rel="stylesheet" href="${request.static_url('magmaweb:static/style.css')}" type="text/css"/>
<style type="text/css">
.x-action-col-cell img.metabolize-col {
	background-image:
		url(${request.extjsroot}/examples/shared/icons/fam/cog.png);
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
	src="${request.static_url('magmaweb:static/d3/d3.min.js')}"></script>
<script type="text/javascript" src="${request.extjsroot}/ext.js"></script>

<script type="text/javascript">
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
      features: {
        // should run buttons (upload ms data, upload molecules, metabolize (one molecule), annotate) be shown
        run: false,
        // should (un)assign buttons be shown
        assign: ${json.dumps(canRun)|n},
        // should logout button be shown
        authenticated: ${json.dumps(request.user is not None)},
        // should login button be hidden and workspace be shown
        anonymous: ${json.dumps(request.registry.settings.get('auto_register', False))|n},
        // should restrictions be applied ie force one spectral tree
        restricted: ${json.dumps(request.registry.settings.get('restricted', False))|n},
      },
      urls: {
        home: '${request.route_path('home')}',
        fragments: '${request.route_path('results',jobid=jobid)}/fragments/{0}/{1}.json',
        mspectra: '${request.route_path('results',jobid=jobid)}/mspectra/{0}.json?mslevel={1}',
        extractedionchromatogram: '${request.route_path('results',jobid=jobid)}/extractedionchromatogram/{0}.json',
        chromatogram: '${request.route_path('results',jobid=jobid)}/chromatogram.json'
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
	<%include file="logos.mak"/>
</body>
</html>
