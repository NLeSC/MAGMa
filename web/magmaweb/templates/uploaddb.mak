<html>
<head>
  <title>MAGMa - Result database upload</title>
  <meta http-equiv="Content-Type" content="text/html;charset=UTF-8"/>
</head>
<body>
<form action="${request.route_url('uploaddb')}" method="post" accept-charset="utf-8"
      enctype="multipart/form-data">

    <label for="db_file">Result database file</label>
    <input id="db_file" name="db_file" type="file" value="" />

    <input type="submit" value="Upload" />
</form></body>
</html>