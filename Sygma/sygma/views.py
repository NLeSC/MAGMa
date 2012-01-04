import uuid, os, json
from sygma.models import DBSession
from sygma.models import Metabolite, Scan, Peak, Fragment, Run, Base
from pyramid.response import Response
from pyramid.view import view_config
from pyramid.httpexceptions import HTTPNotFound, HTTPFound
from sqlalchemy.sql.expression import desc, asc
from sqlalchemy.sql import exists, func
from sqlalchemy import and_
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm import aliased
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
"""Views for pyramid based web application"""

def fetch_job(request):
    """ Fetched job using jobid from request.session[id] or request.params[jobid]"""
    if ('id' in request.session):
        return job_factory(request).fromId(request.session['id'])
    elif ('jobid' in request.params):
        # use request.params['jobid'] to construct job aswell, so job can be bookmarked/shared
        # all following requests use jobid as session['id'], so we dont have to change all urls
        request.session['id'] = request.params['jobid']
        return job_factory(request).fromId(request.params['jobid'])
    else:
        raise HTTPFound(location = request.route_url('home'))

def job_factory(request):
    """ Returns a job factory"""
    from sygma.job import JobFactory
    return JobFactory(request.registry.settings['jobrootdir'], 'results.db')

@view_config(route_name='home', renderer='home.mak')
def home(request):
    """Returns homepage on GET.
     On POST:
     1. creates job session
     2. copies file in 'db' param to job session dir
     3. redirects to results action
    """

    if (request.method == 'POST'):
        # TODO remove results db if it exists

        job = job_factory(request).fromQuery(request.POST['db'].file)
        request.session['id'] = job.id
        return Response(json.dumps({"success": True}), content_type='text/html')

    return dict()

@view_config(route_name='results', renderer='results.mak')
def results(request):
    """Returns results page"""
    job = fetch_job(request)
    return dict(run=job.runInfo(), maxmslevel=job.maxMSLevel())

@view_config(route_name='metabolites.json', renderer='json')
def metabolitesjson(request):
    """Returns json document with metabolites, which can be used in a extjs store

    request.params:

    start
        Offset
    limit
        Maximum nr of metabolites to return
    scanid
        Only return metabolites that have hits in scan with this identifier. Adds score column.
    filter
        Json encoded string which is generated by ExtJS component Ext.ux.grid.FiltersFeature
    sort
        How to sort metabolites. Json encoded string which is an array of objects. Eg.
            [{"property":"probability","direction":"DESC"},{"property":"metid","direction":"ASC"}]
    """
    scanid = request.params['scanid'] if ('scanid' in request.params) else None
    filters = json.loads(request.params['filter']) if ('filter' in request.params) else []
    sorts = json.loads(request.params['sort']) if ('sort' in request.params) else []
    job = fetch_job(request)
    metabolites = job.metabolites(
        start=int(request.params['start']),
        limit=int(request.params['limit']),
        scanid=scanid, filters=filters, sorts=sorts
    )
    scans = job.scansWithMetabolites(scanid=scanid, filters=filters)
    return { 'total':metabolites['total'], 'rows':metabolites['rows'], 'scans':scans}

@view_config(route_name='chromatogram.json', renderer='json')
def chromatogramjson(request):
    """Returns json object with the id, rt and basepeakintensity for each lvl1 scan"""
    job = fetch_job(request)
    return job.chromatogram()

@view_config(route_name='mspectra.json', renderer='json')
def mspectrajson(request):
    """Returns json object with peaks of a scan

    Also returns the cutoff applied to the scan
    and mslevel, precursor.id (parent scan id) and precursor.mz

    request.matchdict['scanid']
        Scan identifier of scan of which to return the mspectra

    request.params.mslevel
        Ms level on which the scan must be. Optional.

    """
    job = fetch_job(request)
    scanid = request.matchdict['scanid']
    mslevel = None
    if ('mslevel' in request.params):
        mslevel = request.params['mslevel']
    from sygma.job import ScanNotFound
    try:
        mspectra = job.mspectra(scanid, mslevel)
        return mspectra
    except ScanNotFound:
        raise HTTPNotFound()

@view_config(route_name='extractedionchromatogram.json', renderer='json')
def extractedionchromatogram(request):
    """Returns json object with the extracted ion chromatogram for a metabolite and the id,rt of scans which have metabolite hits

    request.matchdict['metid']
        Metabolite identifier
    """
    metid = request.matchdict['metid']
    job = fetch_job(request)
    return {
        'chromatogram': job.extractedIonChromatogram(metid),
        'scans': job.scansWithMetabolites(metid=metid )
    }

@view_config(route_name='fragments.json', renderer='json')
def fragments(request):
    """Returns json object with metabolites and its lvl2 fragments when node is not set
    When node is set then returns the children fragments which have node as parent fragment

    Can be used in a Extjs.data.TreeStore.

    request.matchdict['scanid']
        Fragments on scan with this identifier

    request.matchdict['metid']
        Fragments of metabolite with this identifier

    request.params['node']
        The fragment identifier to fetch children fragments for.
    """
    job = fetch_job(request)
    from sygma.job import FragmentNotFound
    try:
        fragments = job.fragments(
            scanid=request.matchdict['scanid'],
            metid=request.matchdict['metid'],
            node=request.params['node']
        )
        return fragments
    except FragmentNotFound:
        raise HTTPNotFound()
