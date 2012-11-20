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
        <input name="displayname" value="${request.user.displayname}" size=50></input>
        </div><div>
        <label for="email">Email</label>
        <input name="email" value="${request.user.email}" size=50></input>
        </div>
        <button>Update</button>
    </form>
    <h1>Jobs</h1>
    <table>
<thead><tr><th>Description</th><th>MS filename</th><th>Created at</th></tr></thead>
<tbody>
        % for job in jobs:
        <tr><td><a href="${request.route_url('results',jobid=job['id'])}">${job['description']}</a></td>
        <td>${job['ms_filename']}</td><td>${job['created_at']}</td>
        </tr>
        % endfor
    </tbody></table>
    <hr></hr>
    <a href="${request.route_url('home')}">Home</a>
</body>
</html>