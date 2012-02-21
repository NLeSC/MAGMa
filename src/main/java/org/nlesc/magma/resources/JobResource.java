package org.nlesc.magma.resources;

import java.io.IOException;
import java.net.URISyntaxException;

import javax.ws.rs.POST;
import javax.ws.rs.Path;

import org.nlesc.magma.BrokerFactory;
import org.nlesc.magma.JobDescriptionFactory;
import org.nlesc.magma.JobStateListener;
import org.nlesc.magma.entities.JobSubmitRequest;
import org.nlesc.magma.entities.JobSubmitResponse;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATInvocationException;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.resources.Job;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.ResourceBroker;

/**
 * Java class which will be hosted at /job
 *
 * @author Stefan Verhoeven <s.verhoeven@esciencecenter.nl>
 *
 */
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

    /**
     * A POST to /job with submit a job.
     *
     * Example json request:
     * <pre>
     * {
     *     "executable": "/bin/sh",
     *     "arguments": [
     *         "magma.sh",
     *         "allinone",
     *         "-n", "1",
     *         "data.mzxml", "smiles.txt", "results.db"
     *     ],
     *     "jobdir": "/tmp/jobdir/",
     *     "prestaged": [
     *         "/home/stefanv/workspace/magmajobmanager/magma.sh",
     *         "/home/stefanv/workspace/magmajobmanager/Magma-1.1.tar.gz",
     *         "data.mzxml", "smiles.txt"
     *     ],
     *     "poststaged": ["results.db"],
     *     "stderr": "stderr.txt",
     *     "stdout": "stdout.txt"
     * }
     * </pre>
     *
     * @param jobsubmission
     * @return JobSubmitResponse which contains a job identifier
     * @throws GATObjectCreationException
     * @throws IOException
     * @throws URISyntaxException
     * @throws GATInvocationException
     */
    @POST
    public JobSubmitResponse submitJob(JobSubmitRequest jobsubmission)
            throws GATObjectCreationException, IOException, URISyntaxException,
            GATInvocationException {

        JobDescription jd = this.jobdescriptionfactory
                .getJobDescription(jobsubmission);

        JobStateListener cb = new JobStateListener(
                GAT.createFile(jobsubmission.jobdir));
        ResourceBroker broker = this.brokerfactory.getBroker();

        // TODO pre staging files are uploaded during submitjob and takes a long
        // time
        // from home it takes 90s to upload Magma-1.1.tar.gz (26M) and
        // data.mzxml.bz2 (18M)
        // from nlesc it takes 15s, this is acceptable when there is only one
        // jobsubmission running at the same time
        // Possible solutions:
        // - Try to fetch Magma tarball from lfn/srm instead from client
        // - Perform submitjob in a new thread, so response can be returned
        // directly, write 'PRE_STAGING' to job state file
        Job job = broker.submitJob(jd, cb, "job.status");

        // TODO Store somevar[jobid] = state,
        // update it in callback so GET /job/{jobid} returns state
        return new JobSubmitResponse(Integer.toString(job.getJobID()));
    }
}
