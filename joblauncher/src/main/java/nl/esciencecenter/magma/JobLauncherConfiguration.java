package nl.esciencecenter.magma;

import javax.validation.Valid;

import org.hibernate.validator.constraints.NotEmpty;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.yammer.dropwizard.config.Configuration;

public class JobLauncherConfiguration extends Configuration {
	@NotEmpty
	@JsonProperty
	private String mac_algorithm;

	@NotEmpty
	@JsonProperty
	private String mac_id;

	@NotEmpty
	@JsonProperty
	private String mac_key;

	@Valid
	@NotEmpty
	@JsonProperty("gat")
	private GATConfiguration gatConfiguration = new GATConfiguration();

	public GATConfiguration getGatConfiguration() {
		return gatConfiguration;
	}

	public String getMac_algorithm() {
		return mac_algorithm;
	}

	public String getMac_id() {
		return mac_id;
	}

	public String getMac_key() {
		return mac_key;
	}
}
