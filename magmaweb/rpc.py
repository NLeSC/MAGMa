from pyramid.view import view_config, view_defaults
from magmaweb.job import make_job_factory

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

    @view_config(route_name='rpc.add_structures', renderer='jsonhtml')
    def add_structures(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().add_structures(self.request.POST, job.maxMSLevel() > 0)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.add_ms_data', renderer='jsonhtml')
    def add_ms_data(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().add_ms_data(self.request.POST, job.metabolitesTotalCount() > 0)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.metabolize', renderer='json')
    def metabolize(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().metabolize(self.request.POST, job.maxMSLevel() > 0)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.metabolize_one', renderer='json')
    def metabolize_one(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().metabolize_one(self.request.POST, job.maxMSLevel() > 0)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.annotate', renderer='json')
    def annotate(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().annotate(self.request.POST)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return { 'success': True, 'jobid': str(job.id) }

    @view_config(route_name='rpc.allinone', renderer='jsonhtml')
    @view_config(route_name='home', renderer='jsonhtml', request_method='POST')
    def allinone(self):
        job = self.new_job()
        """ stage input files and write job script """
        jobquery = job.jobquery().allinone(self.request.POST)
        """ Submit job script to job manager """
        self.job_factory.submitQuery(jobquery)
        return { 'success': True, 'jobid': str(job.id) }
