package org.nlesc.magma.resources;

import static org.junit.Assert.*;

import java.net.URISyntaxException;

import junit.framework.TestCase;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATInvocationException;
import org.gridlab.gat.GATObjectCreationException;
import org.junit.Test;
import org.nlesc.magma.entities.JobSubmitRequest;
import org.nlesc.magma.entities.JobSubmitResponse;

public class JobResourceTest extends TestCase {
	JobResource resource;

	@Override
	protected void setUp() throws Exception {
		super.setUp();
		System.setProperty("gat.adaptor.path", "/home/stefanv/ibisworkspace/javagat_sara/lib/adaptors");
		resource = new JobResource();
	}

	@Override
	protected void tearDown() throws Exception {
		GAT.end();
		super.tearDown();
	}

	@Test
	public void testSubmitJob() {
		try {
			JobSubmitRequest in = new JobSubmitRequest("/somepath", "mzxml");
			JobSubmitResponse out;
			out = resource.submitJob(in);
			assertEquals("12345", out.jobid);
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
		}

	}

}
