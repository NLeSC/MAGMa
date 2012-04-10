import json
from pyramid.view import view_config, view_defaults
from pyramid.response import Response
from magmaweb.job import make_job_factory

# TODO Do not use _html_response but,
# create renderer derived from json renderer with content_type='text/html'
@view_defaults(renderer='json')
class RpcViews(object):
    # TODO RpcViews is a wrapper around JobQuery, make self smaller
    def __init__(self, request):
        self.request = request
        self.job_factory = make_job_factory(request.registry.settings)

    def job(self):
        """ Job of current request """
        return self.job_factory.fromId(self.request.matchdict['jobid'])

    def new_job(self):
        """ Returns clone of job of current request
         or empty job if current request has no job"""
        try:
            return self.job_factory.cloneJob(self.job())
        except KeyError:
            return self.job_factory.fromScratch()

    def _html_response(self, id):
        return Response(json.dumps({"success": True, "jobid": str(id) }), content_type='text/html')

    @view_config(route_name='rpc.add_structures')
    def add_structures(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().add_structures(self.request.POST, job.maxMSLevel() > 0)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return self._html_response(job.id)

    @view_config(route_name='rpc.add_ms_data')
    def add_ms_data(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().add_ms_data(self.request.POST, job.metabolitesTotalCount() > 0)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return self._html_response(job.id)

    @view_config(route_name='rpc.metabolize')
    def metabolize(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().metabolize(self.request.POST, job.maxMSLevel() > 0)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.metabolize_one')
    def metabolize_one(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().metabolize_one(self.request.POST, job.maxMSLevel() > 0)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.annotate')
    def annotate(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().annotate(self.request.POST)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.allinone')
    def allinone(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().allinone(self.request.POST)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return self._html_response(job.id)
