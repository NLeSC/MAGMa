<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Job status</title>
 % if status == 'STOPPED':
<meta http-equiv="refresh" content="0;url=${request.route_url('results',jobid=jobid)}" />
% else:
<meta http-equiv="refresh" content="2" />
% endif
</head>
<body><h2>Status: ${status}
</h2></body>
</html>
