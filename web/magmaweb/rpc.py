"""Module with remote procedure calls views for the magma web application"""
import json
from pyramid.view import view_config, view_defaults
from pyramid.httpexceptions import HTTPInternalServerError
from colander import Invalid
from colander import SchemaNode
from colander import String
from magmaweb.job import make_job_factory
from magmaweb.job import Job
from magmaweb.job import JobSubmissionError


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

    def submit_query(self, query, job):
        """Submit query to job factory

        Raises a HTTPInternalServerError exception when
        job factory fails to communicate with job manager
        """
        try:
            self.job_factory.submitQuery(query, job)
        except JobSubmissionError:
            body = {'success': False, 'msg': 'Unable to submit query'}
            raise HTTPInternalServerError(body=json.dumps(body))

    def _status_url(self, job):
        return self.request.route_url('status.json', jobid=job.id)

    @view_config(route_name='rpc.add_structures', renderer='jsonhtml')
    def add_structures(self):
        """Create a copy of current job or new job
        and submit an add structures job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery(self._status_url(job))
        has_scans = job.db.maxMSLevel() > 0
        params = self.request.POST
        try:
            jobquery = jobquery.add_structures(params, has_scans)
        except Invalid as e:
            # no structures given
            if has_scans and 'structure_database' in params and params['structure_database']:
                # structures will be added by database lookup during extra annotate
                pass
            else:
                sd = SchemaNode(String(), name='structure_database')
                msg = 'Either structures or structures_file or structure_database must be set'
                e.add(Invalid(sd, msg))
                raise e

        if has_scans and 'structure_database' in params and params['structure_database']:
            # add structure database location when structure_database is selected
            key = 'structure_database.' + params['structure_database']
            str_db_loc = self.request.registry.settings[key]
            jobquery = jobquery.annotate(params, False, str_db_loc)

        self.submit_query(jobquery, job)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.add_ms_data', renderer='jsonhtml')
    def add_ms_data(self):
        """Adds ms data to job

        Create a copy of current job or new job and submit an add ms data job

        Annotation is performed if current job has structures
        """
        job = self.new_job()
        try:
            job.ms_filename = self.request.POST['ms_data_file'].filename
        except AttributeError:
            job.ms_filename = 'Uploaded as text'
        has_metabolites = job.db.metabolitesTotalCount() > 0
        jobquery = job.jobquery(self._status_url(job))
        jobquery = jobquery.add_ms_data(self.request.POST, has_metabolites)
        self.submit_query(jobquery, job)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.metabolize', renderer='json')
    def metabolize(self):
        """Create a copy of current job or new job and submit a metabolize job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery(self._status_url(job))
        has_scans = job.db.maxMSLevel() > 0
        jobquery = jobquery.metabolize(self.request.POST, has_scans)
        self.submit_query(jobquery, job)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.metabolize_one', renderer='json')
    def metabolize_one(self):
        """Metabolizes one structure.

        Create a copy of current job or new job and submit a metabolize one job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        has_scans = job.db.maxMSLevel() > 0
        jobquery = job.jobquery(self._status_url(job))
        jobquery = jobquery.metabolize_one(self.request.POST, has_scans)
        self.submit_query(jobquery, job)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.annotate', renderer='json')
    def annotate(self):
        """Create a copy of current job or new job and submit a annotate job

        """
        job = self.new_job()
        jobquery = job.jobquery(self._status_url(job)).annotate(self.request.POST)
        self.submit_query(jobquery, job)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.set_description', renderer='json')
    def set_description(self):
        """Sets description of current job

        ``description`` of post parameter.
        """
        job = self.job
        job.description = self.request.POST['description']
        return {'success': True, 'jobid': str(job.id)}

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
