<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
<title>MAGMa - User</title>
</head>
<body>
	User settings:
	<form>
		<label for="userid">User id</label>
		<input name="userid" disabled value="${request.user.userid}"></input>
		<label for="displayname">Display Name</label>
		<input name="displayname" value="${request.user.displayname}"></input>
		<label for="email">Email</label>
		<input name="email" value="${request.user.email}"></input>
		<button>Update</button>
	</form>
	<hr></hr>
	<a href="${request.route_url('home')}">Home</a>
</body>
</html>