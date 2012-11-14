import json
from pyramid.response import Response
from pyramid.view import forbidden_view_config
from pyramid.view import view_config
from pyramid.view import view_defaults
from pyramid.httpexceptions import HTTPNotFound, HTTPFound
from pyramid.security import has_permission
from pyramid.security import remember
from pyramid.security import forget
from pyramid.security import NO_PERMISSION_REQUIRED
from magmaweb.job import make_job_factory
from magmaweb.job import Job, JobQuery
from magmaweb.user import User


class Views(object):
    def __init__(self, request):
        self.request = request
        self.job_factory = make_job_factory(request.registry.settings)

    @view_config(route_name='home', renderer='home.mak',
                 permission=NO_PERMISSION_REQUIRED)
    def home(self):
        """Returns homepage on GET. """
        return {}

    @view_config(route_name='startjob', renderer='startjob.mak',
                 request_method='GET', permission='view')
    def startjob(self):
        """Returns startjob on GET. """
        return {}

    @view_config(route_name='defaults.json', renderer="json", permission='view')
    def defaults(self):
        """ Returns defaults settings to run a job"""
        return {'success': True,
                'data': JobQuery.defaults()}

    @view_config(route_name='home', renderer='jsonhtml', request_method='POST')
    def allinone(self):
        owner = self.request.user.userid
        job = self.job_factory.fromScratch(owner)
        job.ms_filename = self.request.POST['ms_data_file'].filename
        jobquery = job.jobquery().allinone(self.request.POST)
        self.job_factory.submitQuery(jobquery)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='uploaddb', renderer='uploaddb.mak', permission='view')
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

    @view_config(route_name='workspace',
                 permission='view',
                 renderer='workspace.mak')
    def workspace(self):
        jobs = []
        owner = self.request.user
        for jobmeta in owner.jobs:
            jobs.append({'id': str(jobmeta.jobid),
                         'description': jobmeta.description,
                         'ms_filename': jobmeta.ms_filename,
                         'created_at': str(jobmeta.created_at)})

        return {'jobs': jobs}

    @view_config(route_name='login', renderer='login.mak')
    def login(self):
        login_url = self.request.route_url('login')
        referrer = self.request.url
        if referrer == login_url:
            # never use the login form itself as came_from
            referrer = self.request.route_url('home')
        came_from = self.request.params.get('came_from', referrer)
        userid = ''
        password = ''

        if self.request.method == 'POST':
            userid = self.request.params['userid']
            password = self.request.params['password']
            user = User.by_id(userid)
            if user is not None and user.validate_password(password):
                headers = remember(self.request, userid)
                return HTTPFound(location=came_from,
                                 headers=headers)

        return dict(
            came_from=came_from,
            userid=userid,
            password=password,
            )

    @forbidden_view_config()
    def forbidden(self):
        """Redirect to login if not logged in or give forbidden exception"""
        if self.request.user is None:
            referrer = self.request.url
            if referrer == self.request.route_url('login'):
                # never use the login form itself as came_from
                referrer = self.request.route_url('home')
            query = {'came_from': referrer}
            location = self.request.route_url('login', _query=query)
            return HTTPFound(location=location)
        else:
            return self.request.exception

    @view_config(route_name='logout')
    def logout(self):
        headers = forget(self.request)
        return HTTPFound(location=self.request.route_url('home'),
                         headers=headers)


@view_defaults(context=Job)
class JobViews(object):
    """Views for pyramid based web application with job"""
    def __init__(self, job, request):
        self.request = request
        self.job = job

    @view_config(route_name='status', renderer='status.mak', permission='run')
    @view_config(route_name='status.json', renderer='json', permission='run')
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
        return dict(status=jobstate, jobid=jobid)

    @view_config(route_name='results', renderer='results.mak', permission='view')
    def results(self):
        """Returns results page"""
        db = self.job.db
        # determine if Run buttons should be shown
        canRun = has_permission('run', self.job, self.request)
        return dict(
                    run=db.runInfo(),
                    maxmslevel=db.maxMSLevel(),
                    jobid=self.job.id,
                    canRun=bool(canRun)  # coerce pyramid.security.Allowed|Denied to boolean
                    )

    @view_config(route_name='metabolites.json', renderer='json', permission='view')
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
                     "smiles" : "O=C(O)C1OC(Oc2c(O)cc(CC3OCCC3)cc2)C(O)C(O)C1O",
                     "level" : 1,
                     "isquery" : false,
                     "molformula" : "C17H22O9",
                     "logp" : -0.615300000000001,
                     "mim" : 370.1263823051,
                     "reactionsequence" : "O-glucuronidation_(aromatic_hydroxyl)"
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
                     "reactionsequence" : "O-glucuronidation_(aromatic_hydroxyl)"
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

    @view_config(route_name='metabolites.csv', permission='view')
    def metabolitescsv(self):
        """Same as metabolitesjson(), but returns csv file

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

    @view_config(route_name='metabolites.sdf', permission='view')
    def metabolitessdf(self):
        """Same as metabolitesjson(), but returns sdf file

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

    @view_config(route_name='mspectra.json', renderer='json', permission='view')
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

    @view_config(route_name='extractedionchromatogram.json', renderer='json', permission='view')
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

    @view_config(route_name='fragments.json', renderer='json', permission='view')
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

    @view_config(route_name='runinfo.json', renderer="json", permission='view')
    def runinfojson(self):
        """ Returns settings used for job run or
        if job has not run the default value

        """
        r = self.job.db.runInfo()
        defaults = Views(self.request).defaults()
        if (r is None):
            return defaults
        else:
            runinfo = defaults['data']
            if r.n_reaction_steps:
                runinfo['n_reaction_steps'] = r.n_reaction_steps
            if r.metabolism_types:
                runinfo['metabolism_types'] = r.metabolism_types.split(',')
            if r.ionisation_mode:
                runinfo['ionisation_mode'] = r.ionisation_mode
            if r.skip_fragmentation:
                runinfo['skip_fragmentation'] = r.skip_fragmentation
            if r.ms_intensity_cutoff:
                runinfo['ms_intensity_cutoff'] = r.ms_intensity_cutoff
            if r.msms_intensity_cutoff:
                runinfo['msms_intensity_cutoff'] = r.msms_intensity_cutoff
            if r.mz_precision:
                runinfo['mz_precision'] = r.mz_precision
            if r.use_all_peaks:
                runinfo['use_all_peaks'] = r.use_all_peaks
            if r.abs_peak_cutoff:
                runinfo['abs_peak_cutoff'] = r.abs_peak_cutoff
            if r.rel_peak_cutoff:
                runinfo['rel_peak_cutoff'] = r.rel_peak_cutoff
            if r.max_ms_level:
                runinfo['max_ms_level'] = r.max_ms_level
            if r.precursor_mz_precision:
                runinfo['precursor_mz_precision'] = r.precursor_mz_precision
            if r.max_broken_bonds:
                runinfo['max_broken_bonds'] = r.max_broken_bonds
            return {'success': True, 'data': runinfo}
