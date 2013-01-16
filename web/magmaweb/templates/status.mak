<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Job status</title>
 % if status == 'STOPPED':
<meta http-equiv="refresh"
	content="0;url=${request.route_url('results',jobid=jobid)}" />
% else:
<meta http-equiv="refresh" content="30" />
% endif
</head>
<body>Job is in '${status}' status. This page will reload every
	30 seconds.
</body>
</html>