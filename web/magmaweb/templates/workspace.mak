<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en">
<head>
<title>MAGMa - Jobs</title>
</head>
<body>
    <h1>User settings</h1>
    <form>
    <div>
        <label for="userid">User id</label>
        <input name="userid" disabled value="${request.user.userid}"></input>
        </div><div>
        <label for="displayname">Display Name</label>
        <input name="displayname" value="${request.user.displayname}"></input>
        </div><div>
        <label for="email">Email</label>
        <input name="email" value="${request.user.email}"></input>
        </div>
        <button>Update</button>
    </form>
    <h1>Jobs</h1>
    <ul>
        % for job in jobs:
        <li><a href="${request.route_url('results',jobid=job['id'])}">${job['description']}
                (${job['id']})</a></li>
        % endfor
    </ul>
    <hr></hr>
    <a href="${request.route_url('home')}">Home</a>
</body>
</html>