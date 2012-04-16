import json
from pyramid.response import Response
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPNotFound, HTTPFound
from magmaweb.job import make_job_factory, JobNotFound

class Views(object):
    """Views for pyramid based web application"""
    def __init__(self, request):
        self.request = request
        self.job_factory = make_job_factory(request.registry.settings)

    def jobid(self):
        """ Returns job identifier of current request from request.matchdict['jobid']"""
        try:
            return self.request.matchdict['jobid']
        except KeyError:
            raise HTTPFound(location = self.request.route_url('home'))

    def job(self):
        """ Fetched job using jobid"""
        try:
            return self.job_factory.fromId(self.jobid())
        except JobNotFound:
            raise HTTPNotFound()

    @view_config(route_name='home', renderer='home.mak', request_method='GET')
    def home(self):
        """Returns homepage on GET. """
        return {}

    @view_config(route_name='uploaddb', renderer='uploaddb.mak')
    def uploaddb(self):
        """Upload a sqlitedb as ``db_file`` param in POST request and redirects to job results page"""
        if (self.request.method == 'POST'):
            job = self.job_factory.fromDb(self.request.POST['db_file'].file)
            return HTTPFound(location=self.request.route_url('results', jobid=job.id))
        else:
            return {}

    @view_config(route_name='jobfromscratch')
    def jobfromscratch(self):
        """ Initializes a new job and redirects to its results page"""
        job = self.job_factory.fromScratch()
        return HTTPFound(location=self.request.route_url('results', jobid=job.id))

    @view_config(route_name='results', renderer='results.mak')
    def results(self):
        """Returns results page"""
        job = self.job()
        return dict(
                    run=job.runInfo(),
                    maxmslevel=job.maxMSLevel(),
                    jobid=job.id
                    )

    @view_config(route_name='status', renderer='status.mak')
    @view_config(route_name='status.json', renderer='json')
    def job_status(self):
        """Returns status of a job

        Example json response:

        .. code-block:: python

            {
                "status" : "RUNNING",
                "jobid" : "b1eee101-dcc6-435e-baa8-d35e688c408e"
            }

        """
        jobstate = self.job_factory.state(self.jobid())
        return dict(status=jobstate, jobid=self.jobid())

    @view_config(route_name='metabolites.json', renderer='json')
    def metabolitesjson(self):
        """Returns json document with metabolites, which can be used in a extjs store

        request.params:

        `start`
            Offset
        `limit`
            Maximum nr of metabolites to return
        `scanid`
            Only return metabolites that have hits in scan with this identifier. Adds score column.
        `filter`
            Json encoded string which is generated by ExtJS component Ext.ux.grid.FiltersFeature
        `sort`
            How to sort metabolites. Json encoded string which is an array of objects. Eg.

        .. code-block:: python

                [{"property":"probability","direction":"DESC"},{"property":"metid","direction":"ASC"}]

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
                     "nhits" : null,
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
                     "nr_scans" : 2,
                     "reactionsequence" : "O-glucuronidation_(aromatic_hydroxyl)"
                  },
                  {
                     "mol" : " molblock ...",
                     "nhits" : null,
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
                     "nr_scans" : 2,
                     "reactionsequence" : "O-glucuronidation_(aromatic_hydroxyl)"
                  }
               ]
            }

        """
        request = self.request
        scanid = request.params['scanid'] if ('scanid' in request.params) else None
        filters = json.loads(request.params['filter']) if ('filter' in request.params) else []
        sorts = json.loads(request.params['sort']) if ('sort' in request.params) else []
        job = self.job()
        metabolites = job.metabolites(
            start=int(request.params['start']),
            limit=int(request.params['limit']),
            scanid=scanid, filters=filters, sorts=sorts
        )
        scans = job.scansWithMetabolites(filters=filters)
        return { 'total':metabolites['total'], 'rows':metabolites['rows'], 'scans':scans}

    @view_config(route_name='metabolites.csv')
    def metabolitescsv(self):
        """ Same as metabolitesjson(), but returns csv file instead of a json document """
        mets = self.metabolitesjson()
        csv = self.job().metabolites2csv(mets['rows'])
        response = Response(content_type='text/csv', body=csv.getvalue())
        # response.app_iter does not work on StringIO, so use response.body
        # response.app_iter = csv
        return response

    @view_config(route_name='chromatogram.json', renderer='json')
    def chromatogramjson(self):
        """Returns json object with the id, rt and basepeakintensity for each lvl1 scan

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
        return self.job().chromatogram()

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
        job = self.job()
        scanid = self.request.matchdict['scanid']
        mslevel = None
        if ('mslevel' in self.request.params):
            mslevel = self.request.params['mslevel']
        from magmaweb.job import ScanNotFound
        try:
            return job.mspectra(scanid, mslevel)
        except ScanNotFound:
            raise HTTPNotFound()

    @view_config(route_name='extractedionchromatogram.json', renderer='json')
    def extractedionchromatogram(self):
        """Returns json object with the extracted ion chromatogram for a metabolite and the id,rt of scans which have metabolite hits

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
        job = self.job()
        return {
            'chromatogram': job.extractedIonChromatogram(metid),
            'scans': job.scansWithMetabolites(metid=metid)
        }

    @view_config(route_name='fragments.json', renderer='json')
    def fragments(self):
        """Returns json object with metabolites and its lvl2 fragments when ``node`` is not set
        When node is set then returns the children fragments which have node as parent fragment

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
                     "atoms" : "0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25",
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
        job = self.job()
        request = self.request
        from magmaweb.job import FragmentNotFound
        try:
            fragments = job.fragments(
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
        response.app_iter = self.job().stderr()
        return response
