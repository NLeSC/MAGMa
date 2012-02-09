package org.nlesc.magma.resources;

import static org.junit.Assert.*;
import junit.framework.TestCase;

import org.junit.Test;
import org.nlesc.magma.Main;
import org.nlesc.magma.entities.JobSubmitRequest;
import org.nlesc.magma.entities.JobSubmitResponse;

import com.sun.jersey.api.client.Client;

public class JobResourceTest extends TestCase {
	JobResource resource;

	@Override
	protected void setUp() throws Exception {
		super.setUp();
		resource = new JobResource();
	}

	@Test
	public void testSubmitJob() {
		JobSubmitRequest in = new JobSubmitRequest("/somepath", "mzxml");
		JobSubmitResponse out = resource.submitJob(in);

		assertEquals("12345", out.jobid);
	}

}
