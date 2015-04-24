<!DOCTYPE html>
<html lang="en">

<head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>MAGMa - Login</title>
    <!-- Latest compiled and minified CSS -->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.4/css/bootstrap.min.css">
</head>

<body>
    <div class="container">
        <div class="row">
            <div class="span12">
                <h1>Login</h1>
                <form action="${request.route_url('login')}" method="post">
                    <label for="userid">Username</label>
                    <input name="userid" value="${userid}"></input>
                    <label for="password">Password</label>
                    <input type="password" name="password" value="${password}"></input>
                    <input type="hidden" name="came_from" value="${came_from}" />
                    <button class="btn" type="submit" name="submit">Log in</button>
                </form>
                <p>
                    To login an account is required.
                </p>
                <p>
                    An account can be requested from <a href="mailto:Lars.Ridder@wur.nl">Lars.Ridder@wur.nl</a>
                </p>
                <h3>Account benefits</h3>
                <p>
                    <ul>
                        <li>
                            Allow annotation of all precursor ions in a mzxml-file by candidate molecules from a database.
                            Anonymous access is limited to a single (level 1) precursor ion.
                        </li>
                        <li>
                            No time limit for calculations.
                        </li>
                        <li>
                            Allow candidate molecules larger than 64 atoms from Pubchem, Kegg or Human metabolome database.
                        </li>
                        <li>
                            Option to assign molecules to level 1 precursor ions.
                        </li>
                    </ul>

                </p>
            </div>
        </div>
    </div>
    <footer class="footer">
        <div class="container">
            <a href="${request.route_url('home')}">Home</a></div>
    </footer>
</body>

</html>
