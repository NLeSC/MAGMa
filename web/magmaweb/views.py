"""Module with views for the magma web application"""
import json
import time
from pyramid.response import Response
from pyramid.view import forbidden_view_config
from pyramid.view import view_config
from pyramid.view import view_defaults
from pyramid.httpexceptions import HTTPNotFound
from pyramid.httpexceptions import HTTPFound
from pyramid.httpexceptions import HTTPInternalServerError
from pyramid.security import has_permission
from pyramid.security import remember
from pyramid.security import forget
from pyramid.security import NO_PERMISSION_REQUIRED
from pyramid.interfaces import IAuthenticationPolicy
from pyramid_macauth import MACAuthenticationPolicy
from colander import Invalid
from magmaweb.job import make_job_factory
from magmaweb.job import Job
from magmaweb.job import JobQuery
from magmaweb.job import JobSubmissionError
from magmaweb.user import User


class Views(object):
    def __init__(self, request):
        self.request = request
        self.job_factory = make_job_factory(request.registry.settings)

    @view_config(route_name='home', renderer='home.mak')
    def home(self):
        """Returns homepage on GET. """
        return {}

    @view_config(route_name='startjob', renderer='startjob.mak',
                 request_method='GET', permission='view')
    def startjob(self):
        """Returns startjob on GET. """
        return {}

    @view_config(route_name='defaults.json', renderer="json",
                 permission='view')
    def defaults(self):
        """ Returns defaults settings to run a job"""
        selection = self.request.params.get('selection')
        return {'success': True,
                'data': JobQuery.defaults(selection)}

    @view_config(route_name='startjob', renderer='jsonhtml',
                 request_method='POST', permission='view')
    def allinone(self):
        owner = self.request.user.userid
        job = self.job_factory.fromScratch(owner)
        try:
            job.ms_filename = self.request.POST['ms_data_file'].filename
        except AttributeError:
            job.ms_filename = 'Uploaded as text'
        import logging
        logger = logging.getLogger('magmaweb')
        logger.info(self.request.POST)
        status_url = self.request.route_url('status.json', jobid=job.id)
        jobquery = job.jobquery(status_url)

        if ('structure_database' in self.request.POST
                and self.request.POST['structure_database']):
            # add structure database location when
            # structure_database is selected
            key = 'structure_database.'
            key += self.request.POST['structure_database']
            str_db_loc = self.request.registry.settings[key]
            jobquery = jobquery.allinone(self.request.POST, str_db_loc)
        else:
            jobquery = jobquery.allinone(self.request.POST)

        try:
            self.job_factory.submitQuery(jobquery, job)
        except JobSubmissionError:
            body = {'success': False, 'msg': 'Unable to submit query'}
            raise HTTPInternalServerError(body=json.dumps(body))

        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='uploaddb', renderer='uploaddb.mak',
                 permission='view')
    def uploaddb(self):
        """Upload a sqlitedb as ``db_file`` param in POST request
        and redirects to job results page"""
        if (self.request.method == 'POST'):
            dbfile = self.request.POST['db_file'].file
            owner = self.request.user.userid
            job = self.job_factory.fromDb(dbfile, owner)
            results = self.request.route_url('results', jobid=job.id)
            return HTTPFound(location=results)
        else:
            return {}

    @view_config(route_name='jobfromscratch', permission='view')
    def jobfromscratch(self):
        """ Initializes a new job and redirects to its results page"""
        owner = self.request.user.userid
        job = self.job_factory.fromScratch(owner)
        results = self.request.route_url('results', jobid=job.id)
        return HTTPFound(location=results)

    @view_config(context=Invalid)
    def failed_validation(self):
        """Catches colander.Invalid exceptions
        and returns a ExtJS form submission response
        """
        body = {'success': False, 'errors': self.request.exception.asdict()}
        return HTTPInternalServerError(body=json.dumps(body))

    @view_config(route_name='workspace',
                 permission='view',
                 renderer='workspace.mak')
    def workspace(self):
        jobs = []
        owner = self.request.user
        for jobmeta in owner.jobs:
            # Force ISO 8601 format without microseconds
            # so web browsers can parse it
            created_at = jobmeta.created_at.replace(microsecond=0).isoformat()
            jobs.append({'id': str(jobmeta.jobid),
                         'url': self.request.route_url('results',
                                                       jobid=jobmeta.jobid),
                         'description': jobmeta.description,
                         'ms_filename': jobmeta.ms_filename,
                         'created_at': created_at,
                         })

        return {'jobs': jobs}

    @view_config(route_name='access_token', renderer='json', permission='view')
    def access_token(self):
        """Issue a MAC-Type Access Token"""
        # Get a reference to the MACAuthenticationPolicy plugin.
        policy = self.request.registry.getUtility(IAuthenticationPolicy)
        policy = policy.get_policy(MACAuthenticationPolicy)

        # Generate a new id and secret key for the current user.
        user = self.request.user.userid
        expires_in = self.request.registry.settings['access_token.expires_in']
        expires = time.time() + float(expires_in)
        mac_id, mac_key = policy.encode_mac_id(self.request, user,
                                               expires=expires)

        self.request.response.cache_control = "no-store"
        return {"acesss_token": mac_id,
                "mac_key": mac_key,
                "expires_in": expires_in,
                "token_type": "mac",
                "mac_algorithm": "hmac-sha-1"}

    @view_config(route_name='login', renderer='login.mak',
                 permission=NO_PERMISSION_REQUIRED)
    @forbidden_view_config(renderer='login.mak')
    def login(self):
        """Login page

        - Authenticated -> forbidden exception
        - Unauthenticated
        -- GET login page + MAC challenge
        -- POST -> authenticated -> redirect came_from
                                    or home if came_from=login
        -- POST -> unauthenticated -> login page

        """
        if self.request.user is not None:
            return self.request.exception
        else:
            login_url = self.request.route_url('login')
            referrer = self.request.url
            if referrer == login_url:
                # never use the login form itself as came_from
                referrer = self.request.route_url('home')
            came_from = self.request.params.get('came_from', referrer)
            userid = ''
            password = ''

            if self.request.method == 'POST':
                userid = self.request.POST['userid']
                password = self.request.POST['password']
                user = User.by_id(userid)
                if user is not None and user.validate_password(password):
                    headers = remember(self.request, userid)
                    return HTTPFound(location=came_from,
                                     headers=headers)
            else:
                self.request.response.status_int = 401
                # Add MAC challenge
                self.request.response.headers["WWW-Authenticate"] = "MAC"

            return dict(came_from=came_from,
                        userid=userid,
                        password=password,
                        )

    @view_config(route_name='logout', permission='view')
    def logout(self):
        headers = forget(self.request)
        home = self.request.route_url('home')
        return HTTPFound(location=home, headers=headers)


@view_defaults(context=Job, permission='view')
class JobViews(object):
    """Views for pyramid based web application with job"""
    def __init__(self, job, request):
        self.request = request
        self.job = job

    @view_config(route_name='status',
                 renderer='status.mak',
                 permission='run')
    @view_config(route_name='status.json',
                 renderer='json',
                 permission='run',
                 request_method='GET')
    def job_status(self):
        """Returns status of a job

        Example json response:

        .. code-block:: python

            {
                "status" : "RUNNING",
                "jobid" : "b1eee101-dcc6-435e-baa8-d35e688c408e"
            }

        """
        jobid = self.job.id
        jobstate = self.job.state
        return dict(status=jobstate, jobid=str(jobid))

    @view_config(route_name='status.json', renderer='json',
                 permission='monitor',
                 request_method='PUT')
    def set_job_status(self):
        jobid = self.job.id
        jobstate = self.request.body
        self.job.state = jobstate
        return dict(status=jobstate, jobid=str(jobid))

    @view_config(route_name='results',
                 renderer='results.mak',
                 request_method='GET'
                 )
    def results(self):
        """Returns results page"""
        db = self.job.db
        # determine if Run buttons should be shown
        canRun = has_permission('run', self.job, self.request)
        return dict(run=db.runInfo(),
                    maxmslevel=db.maxMSLevel(),
                    jobid=self.job.id,
                    # coerce pyramid.security.Allowed|Denied to boolean
                    canRun=bool(canRun),
                    job=self.job,
                    )

    @view_config(route_name='results',
                 renderer='json',
                 request_method='PUT',
                 permission='run',
                 )
    def updatejson(self):
        """Update fields of :class:`Job`"""
        body = self.request.json_body
        job = self.job
        job.description = body['description']
        job.ms_filename = body['ms_filename']
        return {'success': True, 'message': 'Updated job'}

    @view_config(route_name='results',
                 renderer='json',
                 request_method='DELETE',
                 permission='run',
                 )
    def deletejson(self):
        """Delete job from user database and deletes job directory"""
        self.job.delete()
        del self.job
        self.request.response.status_int = 204
        return {'success': True, 'message': 'Deleted job'}

    @view_config(route_name='metabolites.json', renderer='json')
    def metabolitesjson(self):
        """Returns json document with metabolites,
        which can be used in a extjs store

        request.params:

        `start`
            Offset
        `limit`
            Maximum nr of metabolites to return
        `scanid`
            Only return metabolites that have hits in scan with this identifier
            Adds score column.
        `filter`
            Json encoded string which is generated by
            ExtJS component Ext.ux.grid.FiltersFeature
        `sort`
            How to sort metabolites.
            Json encoded string which is an array of objects. Eg.

        .. code-block:: python

                [{"property":"probability","direction":"DESC"},
                 {"property":"metid","direction":"ASC"}]

        Example response:

        .. code-block:: python

            {
               "scans" : [
                  {
                     "id" : 1787,
                     "rt" : 42.6626666666667
                  },
                  {
                     "id" : 1789,
                     "rt" : 42.7061666666667
                  }
               ],
               "total" : 2,
               "rows" : [
                  {
                     "mol" : "molblock ...",
                     "nhits" : 2,
                     "metid" : 23,
                     "probability" : 0.248155,
                     "origin" : "5-(3,4)-dihydroxyphenyl-g-valerolactone (F)",
                     "score" : 3,
                     "smiles" : "O=C(O)C1OC(Oc2c(O)cc(CC3OCCC3)...",
                     "level" : 1,
                     "isquery" : false,
                     "molformula" : "C17H22O9",
                     "logp" : -0.615300000000001,
                     "mim" : 370.1263823051,
                     "reactionsequence" : "O-glucuronidation_..."
                  },
                  {
                     "mol" : " molblock ...",
                     "nhits" : 2,
                     "metid" : 24,
                     "probability" : 0.248155,
                     "origin" : "5-(3,4)-dihydroxyphenyl-g-valerolactone (F)",
                     "score" : 3,
                     "smiles" : "O=C(O)C1OC(Oc2cc(CC3OCCC3)ccc2O)C(O)C(O)C1O",
                     "level" : 1,
                     "isquery" : false,
                     "molformula" : "C17H22O9",
                     "logp" : -0.615300000000001,
                     "mim" : 370.1263823051,
                     "reactionsequence" : "O-glucuronidation_..."
                  }
               ]
            }

        """
        request = self.request
        if ('scanid' in request.params):
            scanid = request.params['scanid']
        else:
            scanid = None

        def jd(param):
            return json.loads(request.params[param])

        filters = jd('filter') if ('filter' in request.params) else []
        sorts = jd('sort') if ('sort' in request.params) else []
        job = self.job.db
        metabolites = job.metabolites(
            start=int(request.params['start']),
            limit=int(request.params['limit']),
            scanid=scanid, filters=filters, sorts=sorts
        )
        scans = job.scansWithMetabolites(filters=filters)
        totalUnfiltered = job.metabolitesTotalCount()
        return {'totalUnfiltered': totalUnfiltered,
                'total': metabolites['total'],
                'rows': metabolites['rows'],
                'scans': scans}

    @view_config(route_name='metabolites.csv')
    def metabolitescsv(self):
        """Same as :func:`metabolitesjson`, but returns csv file

        Additional request.params:
          `cols`
            Which metabolite columns should be returned.
            Order of output is same as `cols`.
            A empty list selects all columns.
        """
        mets = self.metabolitesjson()
        cols = []
        if ('cols' in self.request.params):
            cols = json.loads(self.request.params['cols'])
        csv = self.job.db.metabolites2csv(mets['rows'], cols=cols)
        response = Response(content_type='text/csv', body=csv.getvalue())
        # response.app_iter does not work on StringIO, so use response.body
        # response.app_iter = csv
        return response

    @view_config(route_name='metabolites.sdf')
    def metabolitessdf(self):
        """Same as :func:`metabolitesjson`, but returns sdf file

        Additional request.params:
          `cols`
            Which metabolite columns should be returned.
            Order of output is same as `cols`.
            A empty list selects all columns.
        """
        mets = self.metabolitesjson()
        cols = []
        if ('cols' in self.request.params):
            cols = json.loads(self.request.params['cols'])
        sdf = self.job.db.metabolites2sdf(mets['rows'], cols=cols)
        response = Response(content_type='chemical/x-mdl-sdfile',
                            charset='utf-8', body=sdf)
        return response

    @view_config(route_name='chromatogram.json', renderer='json')
    def chromatogramjson(self):
        """Return dict with id, rt and basepeakintensity for each lvl1 scan

        Example response:

        .. code-block:: python

            {
                "scans": [{
                    "rt": 0.013115,
                    "intensity": 14556.6,
                    "id": 1
                }, {
                    "rt": 0.027853333333333334,
                    "intensity": 14144.9,
                    "id": 2
                }],
                "cutoff": 200000.0
            }
        """
        return self.job.db.chromatogram()

    @view_config(route_name='mspectra.json', renderer='json')
    def mspectrajson(self):
        """Returns json object with peaks of a scan

        Also returns the cutoff applied to the scan
        and mslevel, precursor.id (parent scan id) and precursor.mz

        request.matchdict['scanid']
            Scan identifier of scan of which to return the mspectra

        request.params.mslevel
            Ms level on which the scan must be. Optional.

        Example response:

        .. code-block:: python

            {
               "cutoff" : 200000,
               "precursor" : {
                  "mz" : null,
                  "id" : 0
               },
               "mslevel" : 1,
               "peaks" : [
                  {
                     "mz" : 113.024574279785,
                     "intensity" : 32167.611328125
                  },
                  {
                     "mz" : 128.03547668457,
                     "intensity" : 53636.15625
                  }
               ]
            }

        """
        scanid = self.request.matchdict['scanid']
        mslevel = None
        if ('mslevel' in self.request.params):
            mslevel = self.request.params['mslevel']
        from magmaweb.job import ScanNotFound
        try:
            return self.job.db.mspectra(scanid, mslevel)
        except ScanNotFound:
            raise HTTPNotFound()

    @view_config(route_name='extractedionchromatogram.json', renderer='json')
    def extractedionchromatogram(self):
        """Returns json object with the extracted ion chromatogram
        for a metabolite and the id,rt of scans which have metabolite hits

        request.matchdict['metid']
            Metabolite identifier

        Example response:

        .. code-block:: python

             {
               "scans" : [
                  {
                     "id" : 1787,
                     "rt" : 42.6626666666667
                  },
                  {
                     "id" : 1789,
                     "rt" : 42.7061666666667
                  }
               ],
               "chromatogram" : [
                  {
                     "intensity" : 0,
                     "rt" : 0.013115
                  },
                  {
                     "intensity" : 0,
                     "rt" : 0.0278533333333333
                  }
               ]
            }

        """
        metid = self.request.matchdict['metid']
        return {
            'chromatogram': self.job.db.extractedIonChromatogram(metid),
            'scans': self.job.db.scansWithMetabolites(metid=metid)
        }

    @view_config(route_name='fragments.json', renderer='json')
    def fragments(self):
        """Returns json object with metabolites
        and its lvl2 fragments when ``node`` is not set.
        When node is set then returns the children fragments
        which have node as parent fragment.

        Can be used in a Extjs.data.TreeStore.
        From request.matchdict following keys are used:

        ``scanid``
            Fragments on scan with this identifier

        ``metid``
            Fragments of metabolite with this identifier

        ``node``
            The fragment identifier to fetch children fragments for.

        Example response when node='':

        .. code-block:: python

            {
               "expanded" : true,
               "children" : [
                  {
                     "deltah" : -1,
                     "mol" : "molblock ...",
                     "metid" : 23,
                     "fragid" : 5,
                     "score" : 3,
                     "children" : [
                        {
                           "deltah" : -2,
                           "mol" : "molblock ...",
                           "metid" : 23,
                           "fragid" : 6,
                           "score" : 2,
                           "mass" : 115.039519091,
                           "scanid" : 1790,
                           "expanded" : true,
                           "mz" : 113.024360656738,
                           "mslevel" : 2,
                           "atoms" : "14,15,16,20,22,23,24,25",
                           "leaf" : true
                        }
                     ],
                     "mass" : 370.1263823051,
                     "scanid" : 1789,
                     "expanded" : true,
                     "mz" : 369.119262695312,
                     "mslevel" : 1,
                     "atoms" : "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15",
                     "leaf" : false
                  }
               ]
            }

        Example response when node!='':

        .. code-block:: python

                [
                    {
                       "deltah" : -2,
                       "mol" : "molblock ...",
                       "metid" : 23,
                       "fragid" : 6,
                       "score" : 2,
                       "mass" : 115.039519091,
                       "scanid" : 1790,
                       "expanded" : true,
                       "mz" : 113.024360656738,
                       "mslevel" : 2,
                       "atoms" : "14,15,16,20,22,23,24,25",
                       "leaf" : true
                    }
                 ]

        """
        request = self.request
        from magmaweb.job import FragmentNotFound
        try:
            fragments = self.job.db.fragments(
                scanid=request.matchdict['scanid'],
                metid=request.matchdict['metid'],
                node=request.params['node']
            )
            return fragments
        except FragmentNotFound:
            raise HTTPNotFound()

    @view_config(route_name='stderr.txt')
    def stderr(self):
        """Returns file object of stderr.txt file of job"""
        response = Response(content_type='text/plain')
        response.app_iter = self.job.stderr()
        return response

    @view_config(route_name='runinfo.json', renderer="json")
    def runinfojson(self):
        """ Returns settings used for job run or
        if job has not run the default value

        """
        r = self.job.db.runInfo()

        defaults = Views(self.request).defaults()
        if r is None:
            return defaults
        else:
            runinfo = defaults['data']
            for key in r.__dict__.keys():
                if key in runinfo:
                    runinfo[key] = r.__getattribute__(key)
            return {'success': True, 'data': runinfo}
