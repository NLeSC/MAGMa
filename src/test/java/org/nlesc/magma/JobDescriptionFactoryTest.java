package org.nlesc.magma;


import org.gridlab.gat.resources.JobDescription;
import org.junit.Test;
import org.nlesc.magma.entities.JobSubmitRequest;

import junit.framework.TestCase;

public class JobDescriptionFactoryTest extends TestCase {

	@Test
	public void testGetJobDescriptionSleep() {
		JobDescriptionFactory fact = new JobDescriptionFactory();

		JobSubmitRequest jobsubmission = new JobSubmitRequest("/somepath", "sleep");
		try {
			JobDescription jd = fact.getJobDescription(jobsubmission);
			assertEquals(jd.getSoftwareDescription().getExecutable(), "/bin/sleep");
			assertSame(jd.getSoftwareDescription().getArguments()[0], "30");
		} catch (Exception e) {
			fail("Sleep job type must exist");
		}
	}

	@Test
	public void testGetJobDescriptionUnknown() {
		JobDescriptionFactory fact = new JobDescriptionFactory();
		JobSubmitRequest jobsubmission = new JobSubmitRequest("/somepath", "holdbreath");
		try {
			fact.getJobDescription(jobsubmission);
			fail("Holdbreath job type must not exist");
		} catch (Exception e) {
			assertEquals(e.getMessage(), "Unknown job type: 'holdbreath'");
		}
	}

}
