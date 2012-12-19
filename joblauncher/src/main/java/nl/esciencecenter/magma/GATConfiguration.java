package nl.esciencecenter.magma;

import java.util.Set;

import org.gridlab.gat.Preferences;
import org.gridlab.gat.URI;
import org.hibernate.validator.constraints.NotEmpty;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.google.common.collect.ImmutableMap;
import com.yammer.dropwizard.config.Configuration;

public class GATConfiguration {
    @NotEmpty
    @JsonProperty
    private URI brokerURI;

    @JsonProperty
    private ImmutableMap<String, String> preferences = ImmutableMap.of();

    public URI getBrokerURI() {
        return brokerURI;
    }

    public Preferences getPreferences() {
    	Preferences prefs = new Preferences();
    	Set<String> keys = preferences.keySet();
    	for (Object key : keys) {
            prefs.put((String) key, preferences.get(key));
        }
        return prefs;
    }
}
