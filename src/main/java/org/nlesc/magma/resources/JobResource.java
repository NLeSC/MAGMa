package org.nlesc.magma.resources;

import javax.ws.rs.POST;
import javax.ws.rs.Path;

import org.nlesc.magma.BrokerFactory;
import org.nlesc.magma.JobDescriptionFactory;
import org.nlesc.magma.JobStateListener;
import org.nlesc.magma.entities.JobSubmitRequest;
import org.nlesc.magma.entities.JobSubmitResponse;

import org.gridlab.gat.GAT;
import org.gridlab.gat.resources.Job;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.ResourceBroker;

@Path("/job")
public class JobResource {
    protected JobDescriptionFactory jobdescriptionfactory;
    protected BrokerFactory brokerfactory;

    public JobResource() {
        super();
        this.jobdescriptionfactory = new JobDescriptionFactory();
        this.brokerfactory = new BrokerFactory();
    }

    public void setBroker(BrokerFactory broker) {
        this.brokerfactory = broker;
    }

    @POST
    public JobSubmitResponse submitJob(JobSubmitRequest jobsubmission)
            throws Exception {

        JobDescription jd = this.jobdescriptionfactory
                .getJobDescription(jobsubmission);
        ResourceBroker broker = this.brokerfactory.getBroker();

        JobStateListener cb = new JobStateListener(
                GAT.createFile(jobsubmission.jobdir));

        // TODO pre staging files are uploaded during submitjob and takes a long time
        // as it is not done async
        // from home it takes 90s to upload Magma-1.1.tar.gz (26M) and data.mzxml.bz2 (18M)
        // TODO from nlesc ???s
        // Try to fetch Magma tarball from lfn/srm instead from client
        Job job = broker.submitJob(jd, cb, "job.status");

        // TODO Store somevar[jobid] = state,
        // update it in callback so GET /job/{jobid} returns state
        return new JobSubmitResponse(Integer.toString(job.getJobID()));
    }
}
