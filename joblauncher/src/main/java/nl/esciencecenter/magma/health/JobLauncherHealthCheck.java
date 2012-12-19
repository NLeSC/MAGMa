package nl.esciencecenter.magma.health;

import com.yammer.metrics.core.HealthCheck;

public class JobLauncherHealthCheck extends HealthCheck {

	public JobLauncherHealthCheck(String name) {
		super(name);
	}

	@Override
	protected Result check() throws Exception {
		// TODO test if broker is ok
		return Result.healthy();
	}

}
