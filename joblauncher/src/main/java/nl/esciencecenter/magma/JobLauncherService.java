package nl.esciencecenter.magma;

import org.apache.http.client.CredentialsProvider;
import org.apache.http.client.HttpClient;
import org.apache.http.impl.client.AbstractHttpClient;

import nl.esciencecenter.magma.gat.GATManager;
import nl.esciencecenter.magma.health.JobLauncherHealthCheck;
import nl.esciencecenter.magma.mac.MacCredential;
import nl.esciencecenter.magma.mac.MacScheme;
import nl.esciencecenter.magma.mac.MacSchemeFactory;
import nl.esciencecenter.magma.resources.JobLauncherResource;

import com.google.common.collect.ImmutableList;
import com.yammer.dropwizard.Service;
import com.yammer.dropwizard.client.HttpClientBuilder;
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
		HttpClient httpClient = new HttpClientBuilder().using(configuration.getHttpClientConfiguration())
                .build();

		httpClient = MacifyHttpClient((AbstractHttpClient) httpClient, configuration.getMacs());

		environment.addResource(new JobLauncherResource(gatmanager.getBroker(), httpClient));
		environment.addHealthCheck(new JobLauncherHealthCheck("joblauncher"));
	}

	private AbstractHttpClient MacifyHttpClient(AbstractHttpClient httpClient,
			ImmutableList<MacCredential> macs) {

		// Add MAC scheme
		httpClient.getAuthSchemes().register(MacScheme.SCHEME_NAME, new MacSchemeFactory());

		// Add configured MAC id/key pairs.
		CredentialsProvider credentialProvider = httpClient.getCredentialsProvider();
		for (MacCredential mac: macs) {
			credentialProvider.setCredentials(mac.getAuthScope(), mac);
		}

		return httpClient;
	}

}
