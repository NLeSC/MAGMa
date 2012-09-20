<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <title>MAGMa - Job status</title>
  <meta http-equiv="Content-Type" content="text/html;charset=UTF-8"/>
% if status == 'STOPPED':
<meta http-equiv="refresh" content="0;url=${request.route_url('results',jobid=jobid)}"/>
% else:
<meta http-equiv="refresh" content="30"/>
% endif
</head>
<body>
Job is in '${status}' status.
This page will reload every 30 seconds.
</body>
