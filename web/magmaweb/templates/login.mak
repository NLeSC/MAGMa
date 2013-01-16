<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>MAGMa - Login</title>
</head>
<body>
	<h1>Login:</h1>
	<form action="${request.route_url('login')}" method="post">
		<label for="userid">Username</label> <input name="userid"
			value="${userid}"></input> <label for="password">Password</label> <input
			type="password" name="password" value="${password}"></input> <input
			type="hidden" name="came_from" value="${came_from}" />
		<button type="submit" name="submit">Log in</button>
	</form>
	<hr></hr>
	<a href="${request.route_url('home')}">Home</a>
</body>
</html>