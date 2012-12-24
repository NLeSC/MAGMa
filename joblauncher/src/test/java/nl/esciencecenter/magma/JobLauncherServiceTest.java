package nl.esciencecenter.magma;

import static org.junit.Assert.*;

import java.net.URI;
import java.net.URISyntaxException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.http.auth.AuthScope;
import org.apache.http.auth.Credentials;
import org.apache.http.auth.params.AuthPNames;
import org.apache.http.client.params.AuthPolicy;
import org.apache.http.impl.client.DefaultHttpClient;
import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.yammer.dropwizard.client.HttpClientConfiguration;
import com.yammer.dropwizard.config.Bootstrap;
import com.yammer.dropwizard.config.Environment;
import static org.mockito.Matchers.any;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;

import nl.esciencecenter.magma.gat.GATConfiguration;
import nl.esciencecenter.magma.gat.GATManager;
import nl.esciencecenter.magma.health.JobLauncherHealthCheck;
import nl.esciencecenter.magma.mac.MacCredential;
import nl.esciencecenter.magma.mac.MacScheme;
import nl.esciencecenter.magma.resources.JobLauncherResource;

public class JobLauncherServiceTest {
	private final Environment environment = mock(Environment.class);
	private final JobLauncherService service = new JobLauncherService();

	@Test
	public void testInitialize() {
		Bootstrap<JobLauncherConfiguration> bootstrap = new Bootstrap<JobLauncherConfiguration>(
				service);

		service.initialize(bootstrap);

		assertEquals("joblauncher", bootstrap.getName());
	}

	@Test
	public void testRun() throws Exception {
		JobLauncherConfiguration config = sampleConfiguration();

		service.run(config, environment);

		verify(environment).addResource(any(JobLauncherResource.class));
		verify(environment).addHealthCheck(any(JobLauncherHealthCheck.class));
		verify(environment).manage(any(GATManager.class));

		// TODO test injection of MAC Credentials into httpClient
		// or fold injection into extented HttpClientBuilder
	}

	private JobLauncherConfiguration sampleConfiguration()
			throws URISyntaxException {
		ImmutableMap<String, Object> prefs = ImmutableMap.of(
				"localq.max.concurrent.jobs", (Object) 1);
		GATConfiguration gat = new GATConfiguration("localq://localhost", prefs);
		ImmutableList<MacCredential> macs = ImmutableList.of(new MacCredential(
				"id", "key", new URI("http://localhost")));
		HttpClientConfiguration httpClient = new HttpClientConfiguration();

		JobLauncherConfiguration config = new JobLauncherConfiguration(gat,
				macs, httpClient);
		return config;
	}

	@Test
	public void testMacifyHttpClient() throws URISyntaxException {
		JobLauncherConfiguration config = sampleConfiguration();
		DefaultHttpClient httpClient = new DefaultHttpClient();

		JobLauncherService.MacifyHttpClient(httpClient, config.getMacs());

		assertTrue("MAC Registered auth scheme", httpClient.getAuthSchemes()
				.getSchemeNames().contains("mac"));

		MacCredential expected_creds = config.getMacs().get(0);
		AuthScope authscope = expected_creds.getAuthScope();
		Credentials creds = httpClient.getCredentialsProvider().getCredentials(
				authscope);
		assertEquals(expected_creds, creds);

		List<String> authSchemes = Collections
				.unmodifiableList(Arrays.asList(new String[] {
						MacScheme.SCHEME_NAME, AuthPolicy.SPNEGO,
						AuthPolicy.KERBEROS, AuthPolicy.NTLM,
						AuthPolicy.DIGEST, AuthPolicy.BASIC }));
		assertEquals(authSchemes,
				httpClient.getParams()
						.getParameter(AuthPNames.TARGET_AUTH_PREF));
	}
}
