package nl.esciencecenter.magma;

import static com.yammer.dropwizard.testing.JsonHelpers.fromJson;
import static com.yammer.dropwizard.testing.JsonHelpers.jsonFixture;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import nl.esciencecenter.magma.gat.GATConfiguration;
import nl.esciencecenter.magma.mac.MacCredential;

import org.junit.Before;
import org.junit.Test;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.yammer.dropwizard.client.HttpClientConfiguration;

public class JobLauncherConfigurationTest {
	JobLauncherConfiguration conf;

	@Before
	public void testJobLauncherConfiguration() throws URISyntaxException {
		ImmutableMap<String, Object> prefs = ImmutableMap.of("localq.max.concurrent.jobs", (Object) 1);
		GATConfiguration gat = new GATConfiguration("localq://localhost", prefs);
		ImmutableList<MacCredential> macs = ImmutableList.of(new MacCredential("id", "key", new URI("http://localhost")));
		HttpClientConfiguration httpClient = new HttpClientConfiguration();

		conf = new JobLauncherConfiguration(gat, macs, httpClient);

		assertEquals(gat, conf.getGatConfiguration());
		assertEquals(httpClient, conf.getHttpClientConfiguration());
		assertEquals(macs, conf.getMacs());
	}

	@Test
	public void testJobLauncherConfigurationNoArgs() {
		JobLauncherConfiguration conf = new JobLauncherConfiguration();
		assertEquals(new GATConfiguration(), conf.getGatConfiguration());
		assertNotNull(conf.getHttpClientConfiguration());
		ImmutableList<MacCredential> macs = ImmutableList.of();
		assertEquals(macs, conf.getMacs());
	}
	@Test

	public void deserializesFromJSON() throws IOException, URISyntaxException {
		JobLauncherConfiguration conf = fromJson(jsonFixture("fixtures/joblauncher.config.json"),
				JobLauncherConfiguration.class);

		ImmutableMap<String, Object> prefs = ImmutableMap.of("localq.max.concurrent.jobs", (Object) 1);
		GATConfiguration gat = new GATConfiguration("localq://localhost", prefs);
		ImmutableList<MacCredential> macs = ImmutableList.of(new MacCredential("id", "key", new URI("http://localhost")));

		assertEquals(gat, conf.getGatConfiguration());
		assertNotNull(conf.getHttpClientConfiguration());
		assertEquals(macs, conf.getMacs());
	}

}
