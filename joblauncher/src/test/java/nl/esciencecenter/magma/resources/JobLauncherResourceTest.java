package nl.esciencecenter.magma.resources;

import static org.junit.Assert.*;
import static org.mockito.Mockito.*;

import java.net.URI;
import java.net.URISyntaxException;

import nl.esciencecenter.magma.api.JobSubmitRequest;
import nl.esciencecenter.magma.api.JobSubmitResponse;
import nl.esciencecenter.magma.gat.JobStateListener;

import org.apache.http.client.HttpClient;
import org.apache.http.impl.client.DefaultHttpClient;
import org.gridlab.gat.GATInvocationException;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.monitoring.MetricListener;
import org.gridlab.gat.resources.Job;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.ResourceBroker;
import org.junit.Test;

public class JobLauncherResourceTest {

	@Test
	public void testLaunchJob() throws URISyntaxException, GATInvocationException, GATObjectCreationException {
		JobSubmitRequest request = mock(JobSubmitRequest.class);
		URI cb_url = new URI("http://example.com");
		when(request.getStatus_callback_url()).thenReturn(cb_url);
		JobDescription jobdescription = mock(JobDescription.class);
		when(request.toJobDescription()).thenReturn(jobdescription);
		ResourceBroker broker = mock(ResourceBroker.class);
		Job job = mock(Job.class);
		when(job.getJobID()).thenReturn(1234);
		when(broker.submitJob(any(JobDescription.class), any(MetricListener.class), anyString())).thenReturn(job);
		HttpClient httpClient = new DefaultHttpClient();
		JobStateListener listener = new JobStateListener(cb_url, httpClient);

		JobLauncherResource resource = new JobLauncherResource(broker, httpClient);

		JobSubmitResponse response = resource.launchJob(request);

		verify(broker).submitJob(any(JobDescription.class), eq(listener), eq("job.status"));
		assertEquals(response.jobid, "1234");
	}

}
