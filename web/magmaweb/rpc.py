import urllib2
import json
from pyramid.view import view_config, view_defaults
from pyramid.httpexceptions import HTTPInternalServerError
import colander
from magmaweb.job import make_job_factory, Job
from magmaweb.user import status_url


@view_defaults(permission='run', context=Job)
class RpcViews(object):
    """Rpc endpoints"""
    def __init__(self, job, request):
        """View callable with job and request as arguments"""
        self.job = job
        self.request = request
        self.job_factory = make_job_factory(request.registry.settings)

    def new_job(self):
        """Returns clone of job of current request"""
        owner = self.request.user.userid
        job = self.job_factory.cloneJob(self.job, owner)
        return job

    def submit_query(self, query):
        """Submit query to job factory

        Raises a HTTPInternalServerError exception when
        job factory fails to communicate with job manager
        """
        status_callback = status_url(query.id, self.request)
        try:
            self.job_factory.submitQuery(query, status_callback)
        except urllib2.URLError:
            body = {'success': False, 'msg': 'Unable to submit query'}
            raise HTTPInternalServerError(body=json.dumps(body))

    @view_config(route_name='rpc.add_structures', renderer='jsonhtml')
    def add_structures(self):
        """Create a copy of current job or new job
        and submit an add structures job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery().add_structures(self.request.POST,
                                                 job.db.maxMSLevel() > 0)
        self.submit_query(jobquery)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.add_ms_data', renderer='jsonhtml')
    def add_ms_data(self):
        """Adds ms data to job

        Create a copy of current job or new job and submit an add ms data job

        Annotation is performed if current job has structures
        """
        job = self.new_job()
        has_metabolites = job.db.metabolitesTotalCount() > 0
        jobquery = job.jobquery().add_ms_data(self.request.POST,
                                              has_metabolites)
        self.submit_query(jobquery)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.metabolize', renderer='json')
    def metabolize(self):
        """Create a copy of current job or new job and submit a metabolize job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery().metabolize(self.request.POST,
                                             job.db.maxMSLevel() > 0)
        self.submit_query(jobquery)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.metabolize_one', renderer='json')
    def metabolize_one(self):
        """Metabolizes one structure.

        Create a copy of current job or new job and submit a metabolize one job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery().metabolize_one(self.request.POST,
                                                 job.db.maxMSLevel() > 0)
        self.submit_query(jobquery)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.annotate', renderer='json')
    def annotate(self):
        """Create a copy of current job or new job and submit a annotate job

        """
        job = self.new_job()
        jobquery = job.jobquery().annotate(self.request.POST)
        self.submit_query(jobquery)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.set_description', renderer='json')
    def set_description(self):
        """Sets description of current job

        ``description`` of post parameter.
        """
        job = self.job
        job.description = self.request.POST['description']
        return {'success': True, 'jobid': str(job.id)}

    @view_config(context=colander.Invalid)
    def failed_validation(self):
        """Catches colander.Invalid exceptions
        and returns a ExtJS form submission response
        """
        body = {'success': False, 'errors': self.job.asdict()}
        return HTTPInternalServerError(body=json.dumps(body))

    @view_config(route_name='rpc.assign', renderer='json')
    def assign_metabolite2peak(self):
        job = self.job
        scanid = self.request.POST['scanid']
        mz = self.request.POST['mz']
        metid = self.request.POST['metid']
        job.db.assign_metabolite2peak(scanid, mz, metid)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.unassign', renderer='json')
    def unassign_metabolite2peak(self):
        job = self.job
        scanid = self.request.POST['scanid']
        mz = self.request.POST['mz']
        job.db.unassign_metabolite2peak(scanid, mz)
        return {'success': True, 'jobid': str(job.id)}
