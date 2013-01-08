package nl.esciencecenter.magma;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import org.apache.http.auth.params.AuthPNames;
import org.apache.http.client.CredentialsProvider;
import org.apache.http.client.HttpClient;
import org.apache.http.client.params.AuthPolicy;
import org.apache.http.impl.client.AbstractHttpClient;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import nl.esciencecenter.magma.gat.GATManager;
import nl.esciencecenter.magma.gat.JobStateListener;
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

/**
 * Service to submit jobs using a JavaGAT broker.
 *
 * @author verhoes
 *
 */
public class JobLauncherService extends Service<JobLauncherConfiguration> {
    protected static final Logger logger = LoggerFactory
            .getLogger(JobStateListener.class);

    /**
     * Entry point
     *
     * @param args CLI arguments
     * @throws Exception
     */
    public static void main(String[] args) throws Exception {
        new JobLauncherService().run(args);
    }

    @Override
    public void initialize(Bootstrap<JobLauncherConfiguration> bootstrap) {
        bootstrap.setName("joblauncher");
    }

    @Override
    public void run(JobLauncherConfiguration configuration,
            Environment environment) throws Exception {
        GATManager gatmanager = new GATManager(
                configuration.getGatConfiguration());
        environment.manage(gatmanager);
        HttpClient httpClient = new HttpClientBuilder().using(
                configuration.getHttpClientConfiguration()).build();

        httpClient = macifyHttpClient((AbstractHttpClient) httpClient,
                configuration.getMacs());

        environment.addResource(new JobLauncherResource(gatmanager.getBroker(),
                httpClient));
        environment.addHealthCheck(new JobLauncherHealthCheck("joblauncher"));
    }

    /**
     * Adds MAC Access Authentication scheme to http client
     * and registers list of MAC credentials with http client.
     *
     * Http client will use MAC Access Authentication when
     * url is in scope of given MAC credentials.
     *
     * @param httpClient
     * @param macs
     * @return
     */
    public static AbstractHttpClient macifyHttpClient(
            AbstractHttpClient httpClient, ImmutableList<MacCredential> macs) {

        // Add MAC scheme
        httpClient.getAuthSchemes().register(MacScheme.SCHEME_NAME,
                new MacSchemeFactory());

        // Add configured MAC id/key pairs.
        CredentialsProvider credentialProvider = httpClient
                .getCredentialsProvider();
        for (MacCredential mac : macs) {
            credentialProvider.setCredentials(mac.getAuthScope(), mac);
        }

        // Add MAC scheme to ordered list of supported authentication schemes
        // See HTTP authentication parameters chapter on
        // http://hc.apache.org/httpcomponents-client-ga/tutorial/html/authentication.html
        List<String> authSchemes = Collections
                .unmodifiableList(Arrays.asList(new String[] {
                        MacScheme.SCHEME_NAME, AuthPolicy.SPNEGO,
                        AuthPolicy.KERBEROS, AuthPolicy.NTLM,
                        AuthPolicy.DIGEST, AuthPolicy.BASIC }));
        httpClient.getParams().setParameter(AuthPNames.TARGET_AUTH_PREF,
                authSchemes);

        return httpClient;
    }

}
