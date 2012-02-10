package org.nlesc.magma;

import com.sun.jersey.api.container.grizzly2.GrizzlyServerFactory;
import com.sun.jersey.api.core.PackagesResourceConfig;
import com.sun.jersey.api.core.ResourceConfig;
import org.glassfish.grizzly.http.server.HttpServer;

import javax.ws.rs.core.UriBuilder;
import java.io.IOException;
import java.net.URI;

import org.gridlab.gat.GAT;

/**
 * curl -d '{"jobdir":"bla", "jobtype":"foo"}' -H 'Content-Type: application/json' http://localhost:9998/job
 * curl -d '{"jobdir":"/tmp", "jobtype":"mzxmllocal", "arguments":{ "precision":"0.01", "mscutoff":"2e5", "msmscutoff":"0.1", "ionisation":"pos", "nsteps":"2", "phase":"12" }}' -H 'Content-Type: application/json' http://localhost:9998/job
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

	public static void main(String[] args) throws IOException {
	    HttpServer httpServer = startServer();
	    System.out.println(String.format("MaGMA Job manager available at "
	            + "%sapplication.wadl\nTry out %sjob\nHit enter to stop it...",
	            BASE_URI, BASE_URI));
	    System.in.read();
	    httpServer.stop();
		GAT.end();
	}
}
