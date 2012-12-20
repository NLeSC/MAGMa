package nl.esciencecenter.magma.gat;

import org.hibernate.validator.constraints.NotEmpty;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.google.common.collect.ImmutableMap;

public class GATConfiguration {
    @NotEmpty
    @JsonProperty
    private String brokerURI;

    @JsonProperty
    private ImmutableMap<String, Object> preferences = ImmutableMap.of();

	public String getBrokerURI() {
        return brokerURI;
    }

	public void setBrokerURI(String brokerURI) {
		this.brokerURI = brokerURI;
	}

	public ImmutableMap<String, Object> getPreferences() {
		return preferences;
	}

	public void setPreferences(ImmutableMap<String, Object> preferences) {
		this.preferences = preferences;
	}
}
