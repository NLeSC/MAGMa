package org.nlesc.magma.resources;

import static org.mockito.Mockito.*;
import java.net.URISyntaxException;

import junit.framework.TestCase;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATInvocationException;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.monitoring.MetricListener;
import org.gridlab.gat.resources.Job;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.ResourceBroker;
import org.junit.Test;
import org.nlesc.magma.BrokerFactory;
import org.nlesc.magma.JobSubmitCallback;
import org.nlesc.magma.entities.JobSubmitRequest;
import org.nlesc.magma.entities.JobSubmitResponse;

public class JobResourceTest extends TestCase {
	JobResource resource;
	ResourceBroker mockedbroker;

	@Override
	protected void setUp() throws Exception {
		super.setUp();
		resource = new JobResource();
		// mock getBroker() so no job is actually run
		BrokerFactory mockedbrokerfactory = mock(BrokerFactory.class);
		mockedbroker = mock(ResourceBroker.class);
		when(mockedbrokerfactory.getBroker()).thenReturn(mockedbroker);
		resource.setBroker(mockedbrokerfactory);
	}

	@Override
	protected void tearDown() throws Exception {
		GAT.end();
		super.tearDown();
	}

	@Test
	public void testSubmitJob() {
		try {
			Job mockedjob = mock(Job.class);
			when(mockedjob.getJobID()).thenReturn(12345);
			when(
					mockedbroker.submitJob(any(JobDescription.class),
							(MetricListener) anyObject(), anyString()))
					.thenReturn(mockedjob);
			JobSubmitRequest in = new JobSubmitRequest("/somepath", "sleep");

			JobSubmitResponse out = resource.submitJob(in);

			assertEquals("12345", out.jobid);
			// TODO verify jobdescription
			verify(mockedbroker).submitJob(any(JobDescription.class),
					any(JobSubmitCallback.class), eq("job.status"));
		} catch (GATObjectCreationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (GATInvocationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (URISyntaxException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

}
