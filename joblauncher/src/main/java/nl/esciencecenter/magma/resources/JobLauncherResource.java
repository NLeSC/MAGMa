package nl.esciencecenter.magma.resources;

import javax.ws.rs.POST;
import javax.ws.rs.Path;

import org.gridlab.gat.resources.ResourceBroker;

import com.yammer.metrics.annotation.Timed;

import nl.esciencecenter.magma.core.JobSubmitRequest;
import nl.esciencecenter.magma.core.JobSubmitResponse;

@Path("/job")
public class JobLauncherResource {
	protected ResourceBroker broker;

	public JobLauncherResource(ResourceBroker broker) {
		this.broker = broker;
	}

	@POST
	@Timed
	public JobSubmitResponse launchJob(JobSubmitRequest request) {
		return new JobSubmitResponse("submitted");
	}
}
