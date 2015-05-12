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
        self.restricted = self.request.registry.settings['restricted']

    def new_job(self):
        """Returns clone of job of current request"""
        owner = self.request.user.userid
        job = self.job_factory.cloneJob(self.job, owner)
        return job

    def submit_query(self, query, job):
        """Submit query to job factory

        Raises a HTTPInternalServerError exception when
        job factory fails to communicate with job launcher
        """
        try:
            self.job_factory.submitQuery(query, job)
        except JobSubmissionError:
            body = {'success': False, 'msg': 'Unable to submit query'}
            raise HTTPInternalServerError(body=json.dumps(body))

    def _status_url(self, job):
        """Returns status url of `job`"""
        return self.request.route_url('status.json', jobid=job.id)

    @view_config(route_name='rpc.add_structures', renderer='jsonhtml')
    def add_structures(self):
        """Create a copy of current job or new job
        and submit an add structures job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery(self._status_url(job), self.restricted)
        has_scans = job.db.maxMSLevel() > 0
        params = self.request.POST
        try:
            jobquery = jobquery.add_structures(params, has_scans)
        except Invalid as exc:
            # no structures given
            if (has_scans and
                    'structure_database' in params and
                    params['structure_database']):
                # structures will be added by
                # database lookup during extra annotate
                pass
            else:
                node = SchemaNode(String(), name='structure_database')
                msg = 'Either structures or structures_file '
                msg += 'or structure_database must be set'
                exc.add(Invalid(node, msg))
                raise exc

        if (has_scans and
                'structure_database' in params and
                params['structure_database']):
            jobquery = jobquery.annotate(params, False)

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
        has_molecules = job.db.moleculesTotalCount() > 0
        jobquery = job.jobquery(self._status_url(job), self.restricted)
        jobquery = jobquery.add_ms_data(self.request.POST, has_molecules)
        self.submit_query(jobquery, job)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.metabolize', renderer='json')
    def metabolize(self):
        """Create a copy of current job or new job and submit a metabolize job

        Annotation is performed if current job has ms data
        """
        job = self.new_job()
        jobquery = job.jobquery(self._status_url(job), self.restricted)
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
        jobquery = job.jobquery(self._status_url(job), self.restricted)
        jobquery = jobquery.metabolize_one(self.request.POST, has_scans)
        self.submit_query(jobquery, job)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.annotate', renderer='json')
    def annotate(self):
        """Create a copy of current job or new job and submit a annotate job

        """
        job = self.new_job()
        jobquery = job.jobquery(self._status_url(job), self.restricted)
        jobquery = jobquery.annotate(self.request.POST)
        self.submit_query(jobquery, job)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.assign', renderer='json')
    def assign_molecule2peak(self):
        """Assigns molecule with `molid` to peak `mz` in scan `scanid`.
        """
        job = self.job
        scanid = self.request.POST['scanid']
        mz = self.request.POST['mz']
        molid = self.request.POST['molid']
        job.db.assign_molecule2peak(scanid, mz, molid)
        return {'success': True, 'jobid': str(job.id)}

    @view_config(route_name='rpc.unassign', renderer='json')
    def unassign_molecule2peak(self):
        """Unassigns any molecule from peak `mz` in scan `scanid`.
        """
        job = self.job
        scanid = self.request.POST['scanid']
        mz = self.request.POST['mz']
        job.db.unassign_molecule2peak(scanid, mz)
        return {'success': True, 'jobid': str(job.id)}
