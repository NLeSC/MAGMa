package nl.esciencecenter.magma;

import javax.validation.Valid;
import javax.validation.constraints.NotNull;

import nl.esciencecenter.magma.gat.GATConfiguration;
import nl.esciencecenter.magma.mac.MacCredential;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.google.common.collect.ImmutableList;
import com.yammer.dropwizard.client.HttpClientConfiguration;
import com.yammer.dropwizard.config.Configuration;

public class JobLauncherConfiguration extends Configuration {

	@NotNull
	@JsonProperty
	private ImmutableList<MacCredential> macs = ImmutableList.of();

	@Valid
	@NotNull
	@JsonProperty("gat")
	private GATConfiguration gatConfiguration = new GATConfiguration();

	@Valid
	@NotNull
	@JsonProperty
	private HttpClientConfiguration httpClient = new HttpClientConfiguration();

	public JobLauncherConfiguration(GATConfiguration gat,
			ImmutableList<MacCredential> macs,
			HttpClientConfiguration httpClient) {
		this.gatConfiguration = gat;
		this.macs = macs;
		this.httpClient = httpClient;
	}

	public JobLauncherConfiguration() {

	}

	public HttpClientConfiguration getHttpClientConfiguration() {
		return httpClient;
	}

	public GATConfiguration getGatConfiguration() {
		return gatConfiguration;
	}

	public ImmutableList<MacCredential> getMacs() {
		return macs;
	}
}
