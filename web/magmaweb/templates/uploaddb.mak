<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Result database upload</title>
</head>
<body>
	<form action="${request.route_url('uploaddb')}" method="post"
		accept-charset="utf-8" enctype="multipart/form-data">

		<label for="db_file">Result database file</label> <input id="db_file"
			name="db_file" type="file" value="" /> <input type="submit"
			value="Upload" />
	</form>
</body>
</html>