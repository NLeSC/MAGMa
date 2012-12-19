package nl.esciencecenter.magma.resources;

import static org.junit.Assert.*;
import static org.mockito.Mockito.mock;

import java.net.URI;
import java.net.URISyntaxException;

import nl.esciencecenter.magma.core.JobSubmitRequest;
import nl.esciencecenter.magma.core.JobSubmitResponse;

import org.gridlab.gat.resources.ResourceBroker;
import org.junit.Test;

public class JobLauncherResourceTest {

	@Test
	public void testLaunchJob() throws URISyntaxException {
		String[] args = {"magma.sh", "1234"};
		String[] preStaged = {"magma.sh"};
		String[] postStaged = {};
		JobSubmitRequest request = new JobSubmitRequest("/jobdir", "/bin/sh", args,
				preStaged, postStaged,
				"stderr.txt", "stdout.txt",
				new URI("http://localhost/callback")
				);

		ResourceBroker broker = mock(ResourceBroker.class);
		JobLauncherResource resource = new JobLauncherResource(broker);

		JobSubmitResponse response = resource.launchJob(request);

		assertEquals(response.jobid, "submitted");
	}

}
