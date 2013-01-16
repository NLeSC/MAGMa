package nl.esciencecenter.magma;

import javax.validation.Valid;
import javax.validation.constraints.NotNull;

import nl.esciencecenter.magma.gat.GATConfiguration;
import nl.esciencecenter.magma.mac.MacCredential;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.google.common.collect.ImmutableList;
import com.yammer.dropwizard.client.HttpClientConfiguration;
import com.yammer.dropwizard.config.Configuration;

/**
 * Main configuration of job launcher.
 *
 * @author verhoes
 *
 */
public class JobLauncherConfiguration extends Configuration {
    /**
     * List of mac credentials used to perform authentication
     * when using http client on their scope.
     */
    @NotNull
    @JsonProperty
    private ImmutableList<MacCredential> macs = ImmutableList.of();


    @Valid
    @NotNull
    @JsonProperty("gat")
    private GATConfiguration gatConfiguration = new GATConfiguration();

    /**
     * Http client used to PUT state changes of jobs
     * the launcher is monitoring to the status callback url of the job
     */
    @Valid
    @NotNull
    @JsonProperty
    private HttpClientConfiguration httpClient = new HttpClientConfiguration();

    /**
     * Constructor
     *
     * @param gat
     * @param macs
     * @param httpClient
     */
    public JobLauncherConfiguration(GATConfiguration gat,
            ImmutableList<MacCredential> macs,
            HttpClientConfiguration httpClient) {
        this.gatConfiguration = gat;
        this.macs = macs;
        this.httpClient = httpClient;
    }

    /**
     * No argument contructor required for JAXB
     */
    public JobLauncherConfiguration() {

    }

    /**
     *
     * @return HttpClientConfiguration
     */
    public HttpClientConfiguration getHttpClientConfiguration() {
        return httpClient;
    }

    /**
     *
     * @return GATConfiguration
     */
    public GATConfiguration getGatConfiguration() {
        return gatConfiguration;
    }

    /**
     * Credentials used for http client.
     *
     * @return List of MAC Access Authentication credentials
     */
    public ImmutableList<MacCredential> getMacs() {
        return macs;
    }
}
