package nl.esciencecenter.magma.gat;

import org.hibernate.validator.constraints.NotEmpty;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.google.common.base.Objects;
import com.google.common.collect.ImmutableMap;

/**
 * JavaGAT configuration
 *
 * @author verhoes
 *
 */
public class GATConfiguration {
    /**
     * Broker URI used to submit jobs
     */
    @NotEmpty
    @JsonProperty
    private String brokerURI;

    /**
     * JavaGAT preferences, these could also be put javagat.properties file,
     * but I like broker together with it's preferences
     */
    @JsonProperty
    private ImmutableMap<String, Object> preferences = ImmutableMap.of();

    public GATConfiguration(String brokerURI,
            ImmutableMap<String, Object> preferences) {
        this.brokerURI = brokerURI;
        this.preferences = preferences;
    }

    public GATConfiguration() {
        this.brokerURI = null;
    }

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

    @Override
    public int hashCode() {
        return Objects.hashCode(brokerURI, preferences);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        GATConfiguration other = (GATConfiguration) obj;
        return Objects.equal(this.brokerURI, other.brokerURI)
                && Objects.equal(this.preferences, other.preferences);
    }
}