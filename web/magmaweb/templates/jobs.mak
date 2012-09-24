<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <title>MAGMa - Jobs</title>
  </head>
  <body>
  My MAGMa Jobs:
  <ul>
  % for job in jobs:
  <li>
  <a href="${request.route_url('results',jobid=job['id'])}">${job['description']} (${job['id']})</a>
  </li>
  % endfor
  </ul>
  <hr></hr>
  <a href="${request.route_url('home')}">Home</a>
  </body>