package nl.esciencecenter.magma;

import nl.esciencecenter.magma.health.JobLauncherHealthCheck;
import nl.esciencecenter.magma.resources.JobLauncherResource;

import com.yammer.dropwizard.Service;
import com.yammer.dropwizard.config.Bootstrap;
import com.yammer.dropwizard.config.Environment;

public class JobLauncherService extends Service<JobLauncherConfiguration> {
	public static void main(String[] args) throws Exception {
        new JobLauncherService().run(args);
    }

	@Override
	public void initialize(Bootstrap<JobLauncherConfiguration> bootstrap) {
		bootstrap.setName("joblauncher");
	}

	@Override
	public void run(JobLauncherConfiguration configuration, Environment environment)
			throws Exception {
		GATManager gatmanager = new GATManager(configuration.getGatConfiguration());
		environment.manage(gatmanager);
		environment.addResource(new JobLauncherResource(gatmanager.getBroker()));
		environment.addHealthCheck(new JobLauncherHealthCheck("joblauncher"));
	}

}
