package nl.esciencecenter.magma.health;

import static org.junit.Assert.*;

import org.junit.Test;

import com.yammer.metrics.core.HealthCheck.Result;

public class JobLauncherHealthCheckTest {

	@Test
	public void testCheck() throws Exception {
		JobLauncherHealthCheck hc = new JobLauncherHealthCheck("gat");
		assertEquals(Result.healthy(), hc.check());
	}

}
