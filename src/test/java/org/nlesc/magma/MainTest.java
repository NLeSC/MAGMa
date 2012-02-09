package org.nlesc.magma;

import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.MultivaluedMap;

import org.codehaus.jettison.json.JSONException;
import org.codehaus.jettison.json.JSONObject;
import org.glassfish.grizzly.http.server.HttpServer;
import org.nlesc.magma.entities.JobSubmitRequest;

import com.sun.jersey.core.header.MediaTypes;
import com.sun.jersey.core.util.MultivaluedMapImpl;
import com.sun.jersey.api.client.Client;
import com.sun.jersey.api.client.WebResource;

import junit.framework.TestCase;

public class MainTest extends TestCase {

	private HttpServer httpServer;

	private WebResource r;

	public MainTest(String testName) {
		super(testName);
	}

	@Override
	protected void setUp() throws Exception {
		super.setUp();

		// start the Grizzly2 web container
		httpServer = Main.startServer();

		// create the client
		Client c = Client.create();
		r = c.resource(Main.BASE_URI);
	}

	@Override
	protected void tearDown() throws Exception {
		super.tearDown();

		httpServer.stop();
	}

	/**
	 * Test to see that the message "Got it!" is sent in the response.
	 * @throws JSONException
	 */
	public void testJobResource() throws JSONException {
		JSONObject requestMsg = new JSONObject();
		requestMsg.put("jobdir", "/somepath").put("jobtype", "mzxml");

		JSONObject responseMsg = r.path("job").post(JSONObject.class, requestMsg);

		assertEquals("12345", responseMsg.get("jobid"));
	}

	/**
	 * Test if a WADL document is available at the relative path
	 * "application.wadl".
	 */
	public void testApplicationWadl() {
		String serviceWadl = r.path("application.wadl").accept(MediaTypes.WADL)
				.get(String.class);

		assertTrue(serviceWadl.length() > 0);
	}
}
