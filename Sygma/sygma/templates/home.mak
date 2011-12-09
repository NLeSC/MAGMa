<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
  <title>MSygma - Homepage</title>
</head>
<body>
% if 'dbname' in request.session:
<a href="${request.route_url('results')}">Results</a>
% endif
<form action="${request.route_url('home')}" method="post" accept-charset="utf-8"
      enctype="multipart/form-data">

    <label for="db">Sqlite results database</label>
    <input id="db" name="db" type="file" value="" />

    <input type="submit" value="submit" />
</form>
</body>
</html>