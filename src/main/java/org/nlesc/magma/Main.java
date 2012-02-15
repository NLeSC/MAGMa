package org.nlesc.magma;

import com.sun.jersey.api.container.grizzly2.GrizzlyServerFactory;
import com.sun.jersey.api.core.PackagesResourceConfig;
import com.sun.jersey.api.core.ResourceConfig;
import org.glassfish.grizzly.http.server.HttpServer;

import javax.ws.rs.core.UriBuilder;

import java.io.Console;
import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATContext;
import org.gridlab.gat.security.CertificateSecurityContext;

/**
 * REST webservice to submit magma jobs to grid.
 *
 * mkdir /tmp/jobdir
 * put inputfiles in /tmp/jobdir
 * curl -d '{"jobdir":"/tmp/jobdir", "jobtype":"sleep"}' -H 'Content-Type: application/json' http://localhost:9998/job
 * curl -d '{"jobdir":"/tmp/jobdir", "jobtype":"mzxmllocal", "arguments":{ "precision":"0.01", "mscutoff":"2e5", "msmscutoff":"0.1", "ionisation":"1", "nsteps":"2", "phase":"12" }}' -H 'Content-Type: application/json' http://localhost:9998/job
 * curl -d '{"jobdir":"/tmp/jobdir", "jobtype":"mzxmlremote", "arguments":{ "precision":"0.01", "mscutoff":"2e5", "msmscutoff":"0.1", "ionisation":"1", "nsteps":"2", "phase":"12" }}' -H 'Content-Type: application/json' http://localhost:9998/job
 *
 * @author stefanv
 *
 */
public class Main {

    private static URI getBaseURI() {
        // TODO move port to config file
        return UriBuilder.fromUri("http://localhost/").port(9998).build();
    }

    public static final URI BASE_URI = getBaseURI();

    protected static HttpServer startServer() throws IOException {
        System.out.println("Starting grizzly...");
        ResourceConfig rc = new PackagesResourceConfig("org.nlesc.magma.resources");
        return GrizzlyServerFactory.createHttpServer(BASE_URI, rc);
    }

    /**
     * Asks for passphrase and starts webserver.
     *
     * @param args
     * @throws IOException
     * @throws URISyntaxException
     */
    public static void main(String[] args) throws IOException, URISyntaxException {
        setupGATContext();
        HttpServer httpServer = startServer();
        System.out.println(String.format("MaGMA Job manager available at "
                + "%sapplication.wadl%nTry out %sjob%nHit enter to stop it...",
                BASE_URI, BASE_URI));
        System.in.read();
        httpServer.stop();
        GAT.end();
    }

    // Ask the user for the password needed to perform grid-proxy-init
    protected static String getPassphrase() {
        // console() does not work within Eclipse
        Console cons = System.console();
        return new String(cons.readPassword("Enter password for %s (leave empty to use localhost broker):", "Glite certificate"));
    }

    protected static void setupGATContext() throws URISyntaxException {
        CertificateSecurityContext securityContext = new CertificateSecurityContext(
                new org.gridlab.gat.URI(System.getProperty("user.home") + "/.globus/userkey.pem"),
                new org.gridlab.gat.URI(System.getProperty("user.home") + "/.globus/usercert.pem"),
                getPassphrase());

        // Store this SecurityContext in a GATContext
        GATContext context = new GATContext();
        context.addSecurityContext(securityContext);

        context.addPreference("VirtualOrganisation", "nlesc.nl");
        context.addPreference("vomsServerUrl", "voms.grid.sara.nl");
        context.addPreference("vomsServerPort", "30025");
        context.addPreference("vomsHostDN", "/O=dutchgrid/O=hosts/OU=sara.nl/CN=voms.grid.sara.nl");
        context.addPreference("LfcServer", "lfc.grid.sara.nl");
        context.addPreference("bdiiURI", "ldap://bdii.grid.sara.nl:2170");

        GAT.setDefaultGATContext(context);
    }
}
