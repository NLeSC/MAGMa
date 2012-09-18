import urllib2
import json
from pyramid.view import view_config, view_defaults
from pyramid.httpexceptions import HTTPInternalServerError
import colander
from magmaweb.job import make_job_factory

@view_defaults(permission='run')
class RpcViews(object):
    """Rpc endpoints"""
    def __init__(self, job, request):
        """ View callable with request as argument or context, request as arguments"""
        self.job = job
        self.request = request
        self.job_factory = make_job_factory(request.registry.settings)

    def new_job(self):
        """Returns clone of job of current request"""
        return self.job_factory.cloneJob(self.job)

    def submit_query(self, query):
        """ Submit query to job factory

        Raises a HTTPInternalServerError exception when job factory fails to communicate with job manager
        """
        try:
            self.job_factory.submitQuery(query)
        except urllib2.URLError:
            body = { 'success': False, 'msg': 'Unable to submit query'}
            raise HTTPInternalServerError(body=json.dumps(body))

    @view_config(route_name='rpc.add_structures', renderer='jsonhtml')
    def add_structures(self):
        """Create a copy of current job or new job and submit an add structures job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery().add_structures(self.request.POST, job.maxMSLevel() > 0)
        self.submit_query(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.add_ms_data', renderer='jsonhtml')
    def add_ms_data(self):
        """Create a copy of current job or new job and submit an add ms data job

        Annotation is performed if current job has structures
        """
        job = self.new_job()
        jobquery = job.jobquery().add_ms_data(self.request.POST, job.metabolitesTotalCount() > 0)
        self.submit_query(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.metabolize', renderer='json')
    def metabolize(self):
        """Create a copy of current job or new job and submit a metabolize job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery().metabolize(self.request.POST, job.maxMSLevel() > 0)
        self.submit_query(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.metabolize_one', renderer='json')
    def metabolize_one(self):
        """Create a copy of current job or new job and submit a metabolize one job
        Metabolizes one structure.

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery().metabolize_one(self.request.POST, job.maxMSLevel() > 0)
        self.submit_query(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.annotate', renderer='json')
    def annotate(self):
        """Create a copy of current job or new job and submit a annotate job

        """
        job = self.new_job()
        jobquery = job.jobquery().annotate(self.request.POST)
        self.submit_query(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.allinone', renderer='jsonhtml')
    @view_config(route_name='home', renderer='jsonhtml', request_method='POST')
    def allinone(self):
        """Create a copy of current job or new job and submit a allinone job

        Adds structures and ms data, then metabolizes all structures and annotates structures/peaks
        """
        job = self.new_job()
        jobquery = job.jobquery().allinone(self.request.POST)
        self.submit_query(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.set_description', renderer='json')
    def set_description(self):
        """Sets description of current job

        ``description`` of post parameter.
        """
        job = self.job
        desc = self.request.POST['description']
        job.description(desc)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(context=colander.Invalid)
    def failed_validation(self):
        """ Catches colander.Invalid exceptions and returns a ExtJS form submission response"""
        body = {'success': False, 'errors': self.job.asdict()}
        return HTTPInternalServerError(body=json.dumps(body))

    @view_config(route_name='rpc.assign', renderer='json')
    def assign_metabolite2peak(self):
        job = self.job
        scanid = self.request.POST['scanid']
        mz = self.request.POST['mz']
        metid = self.request.POST['metid']
        job.assign_metabolite2peak(scanid, mz, metid)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.unassign', renderer='json')
    def unassign_metabolite2peak(self):
        job = self.job
        scanid = self.request.POST['scanid']
        mz = self.request.POST['mz']
        job.unassign_metabolite2peak(scanid, mz)
        return { 'success': True, 'jobid': str(job.id) }

