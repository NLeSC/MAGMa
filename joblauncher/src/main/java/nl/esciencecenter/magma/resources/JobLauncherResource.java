package nl.esciencecenter.magma.resources;

import javax.ws.rs.POST;
import javax.ws.rs.Path;

import org.apache.http.client.HttpClient;
import org.gridlab.gat.GATInvocationException;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.monitoring.MetricListener;
import org.gridlab.gat.resources.Job;
import org.gridlab.gat.resources.ResourceBroker;

import com.yammer.metrics.annotation.Timed;

import nl.esciencecenter.magma.api.JobSubmitRequest;
import nl.esciencecenter.magma.api.JobSubmitResponse;
import nl.esciencecenter.magma.gat.JobStateListener;

@Path("/job")
public class JobLauncherResource {
	protected ResourceBroker broker;
	protected HttpClient httpClient;

	public JobLauncherResource(ResourceBroker broker, HttpClient httpClient) {
		this.broker = broker;
		this.httpClient = httpClient;
	}

	@POST
	@Timed
	public JobSubmitResponse launchJob(JobSubmitRequest request)
			throws GATInvocationException, GATObjectCreationException {
		MetricListener listener = new JobStateListener(
				request.getStatus_callback_url(), httpClient);

		Job job = broker.submitJob(request.getJobDescription(), listener,
				"job.status");

		return new JobSubmitResponse(Integer.toString(job.getJobID()));
	}
}
