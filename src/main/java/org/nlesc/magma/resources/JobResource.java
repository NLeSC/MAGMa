package org.nlesc.magma.resources;
import javax.ws.rs.Consumes;
import javax.ws.rs.POST;
import javax.ws.rs.Path;
import javax.ws.rs.Produces;

import org.nlesc.magma.entities.JobSubmitRequest;
import org.nlesc.magma.entities.JobSubmitResponse;

@Path("/job")
public class JobResource {
	@POST
	@Produces("application/json")
	@Consumes("application/json")
	public JobSubmitResponse submitJob(JobSubmitRequest jobsubmission) {
		System.err.println(jobsubmission.jobdir + " + " + jobsubmission.jobtype);
		return new JobSubmitResponse("12345");
	}
}
