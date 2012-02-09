package org.nlesc.magma.resources;
import java.net.URISyntaxException;

import javax.ws.rs.POST;
import javax.ws.rs.Path;

import org.nlesc.magma.JobSubmitCallback;
import org.nlesc.magma.entities.JobSubmitRequest;
import org.nlesc.magma.entities.JobSubmitResponse;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATInvocationException;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.URI;
import org.gridlab.gat.resources.Job;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.ResourceBroker;
import org.gridlab.gat.resources.SoftwareDescription;

@Path("/job")
public class JobResource {
	@POST
	public JobSubmitResponse submitJob(JobSubmitRequest jobsubmission) throws GATObjectCreationException, URISyntaxException, GATInvocationException, InterruptedException {

		SoftwareDescription sd = new SoftwareDescription();
		// simulate job by sleeping a while
        sd.setExecutable("/bin/sleep");
        String[] args = { "30" };
        sd.setArguments(args);
//        File stdout = GAT.createFile("hostname.txt");
//        sd.setStdout(stdout);

        JobDescription jd = new JobDescription(sd);
        ResourceBroker broker = GAT.createResourceBroker(new URI(
                "any://localhost"));

        JobSubmitCallback cb = new JobSubmitCallback(jobsubmission.jobdir);
        Job job = broker.submitJob(jd, cb, "job.status");

        // use callback to finalize job
//        while ((job.getState() != JobState.STOPPED)
//                && (job.getState() != JobState.SUBMISSION_ERROR))
//            Thread.sleep(1000);
//
//		System.err.println(jobsubmission.jobdir + " + " + jobsubmission.jobtype + " == " + job.getJobID());

		String jobidstr = Integer.toString(job.getJobID());
		return new JobSubmitResponse(jobidstr);
	}
}
