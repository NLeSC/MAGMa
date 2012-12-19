package nl.esciencecenter.magma;

import org.gridlab.gat.URI;
import org.hibernate.validator.constraints.NotEmpty;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.google.common.collect.ImmutableMap;

public class GATConfiguration {
    @NotEmpty
    @JsonProperty
    private URI brokerURI;

    @JsonProperty
    private ImmutableMap<String, String> preferences = ImmutableMap.of();

    public URI getBrokerURI() {
        return brokerURI;
    }

	public ImmutableMap<String, String> getPreferences() {
		return preferences;
	}

	public void setPreferences(ImmutableMap<String, String> preferences) {
		this.preferences = preferences;
	}
}
