<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
<title>MAGMa - Login</title>
</head>
<body>
	<h1>Login:</h1>
	<form action="${request.route_url('login')}" method="post">
		<label for="userid">Username</label>
		<input name="userid" value="${userid}"></input>
		<label for="password">Password</label>
		<input type="password" name="password" value="${password}"></input>
		<input type="hidden" name="came_from" value="${came_from}" />
		<button type="submit" name="submit">Log in</button>
	</form>
	<hr></hr>
	<a href="${request.route_url('home')}">Home</a>
</body>
</html>