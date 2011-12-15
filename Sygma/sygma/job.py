from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql import func
from sygma.models import Metabolite, Scan, Peak, Fragment, Run, Base
import uuid, os

class JobFactory:
    """ Factory which can create jobs """
    dbname = 'results.db'
    def __init__(self, jobrootdir, dbname):
        """
        jobrootdir
            Directory in which jobs are created, retrieved

        dbname
            Sqlite db file name in job directory
        """
        self.jobrootdir = jobrootdir
        self.dbname = dbname

    def fromId(self, jobid):
        """
        Finds job db in job root dir.
        Returns a Job instance

        jobid
            Job identifier
        """
        session = sessionmaker(bind=create_engine(self.id2url(jobid)))
        return Job(jobid, session())

    def fromQuery(self, dbfile):
        """
        A job directory is created and the dbfile is copied into it
        Returns a Job instance

        dbfile
            The sqlite result db

        TODO replace with arguments with submit form values like mzxml file, metabolite smiles and sygma config
        """
        jobid = uuid.uuid1()

        # create job dir
        os.makedirs(self.id2jobdir(jobid))
        # copy dbfile into job dir
        jobdb = open(self.id2db(jobid), 'wb')
        dbfile.seek(0)
        while 1:
            data = dbfile.read(2<<16)
            if not data:
               break
            jobdb.write(data)
        jobdb.close()

        return self.fromId(jobid)

    def id2jobdir(self, id):
        """ Returns job directory based on id and jobrootdir """
        return os.path.join(self.jobrootdir, str(id))

    def id2db(self, id):
        """ Returns sqlite db of job with id """
        return os.path.join(self.id2jobdir(id), self.dbname)

    def id2url(self, id):
        """ Returns sqlalchemy url of sqlite db of job with id """
        # 3rd / is for username:pw@host which sqlite does not need
        return 'sqlite:///'+self.id2db(id)

class Job:
    """
    Job contains query and results of MSygma calculation run
    """
    def __init__(self, id, session):
        """
        id
            A UUID of the job
        session
            Sqlalchemy session to read/write to job db
        """
        self.id = id
        self.session = session

    def maxMSLevel(self):
        """ Returns the maximum nr of MS levels """
        return self.session.query(func.max(Scan.mslevel)).scalar()

    def runInfo(self):
        """ Returns run info"""
        return self.session.query(Run).one()

