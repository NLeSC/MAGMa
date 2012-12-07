package nl.nlesc.magma;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.net.URI;
import java.net.URISyntaxException;
import java.util.UUID;

import org.codehaus.jettison.json.JSONException;
import org.codehaus.jettison.json.JSONObject;
import org.gridlab.gat.GAT;
import org.gridlab.gat.GATInvocationException;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.io.File;
import org.junit.Test;

import com.sun.jersey.core.header.MediaTypes;
import com.sun.jersey.test.framework.JerseyTest;
import com.sun.jersey.api.client.WebResource;

public class MainTest extends JerseyTest {
	public MainTest() throws Exception {
		super("nl.nlesc.magma.resources");
	}

	/**
	 * Test to see if a job can be submitted This is a integration test, testing
	 * webservice and javagat localq job broker together
	 *
	 * @throws JSONException
	 * @throws GATObjectCreationException
	 * @throws GATInvocationException
	 * @throws URISyntaxException
	 * @throws InterruptedException
	 */
	@Test
	public void testJobResourcePOST() throws JSONException,
			GATObjectCreationException, GATInvocationException,
			URISyntaxException, InterruptedException {
		String jobid = UUID.randomUUID().toString();
		String jobdir = System.getProperty("java.io.tmpdir") + "/" + jobid;
		File jobdirfile = GAT.createFile(jobdir);
		jobdirfile.mkdir();
		// use tests own web server as status callback url
		URI status_cb_url = getBaseURI();

		JSONObject requestMsg = new JSONObject();
		String[] args = { "1" };
		requestMsg.put("jobdir", jobdir);
		requestMsg.put("executable", "sleep");
		requestMsg.put("arguments", args);
		requestMsg.put("status_callback_url", status_cb_url);
//		requestMsg = new JSONObject("{\"jobdir\" : \"/tmp/myjob/\", \"status_callback_url\" : \"http://localhost/job/myjob/status\", \"executable\" : \"/usr/bin/uptime\"}");

		WebResource webResource = resource();
		JSONObject responseMsg = webResource.path("job").post(JSONObject.class,
				requestMsg);



		assertEquals("0", responseMsg.get("jobid")); // first job == 0 and GAT
														// is restarted each
														// time

		// let job finish
		Thread.sleep(2000);

		// TODO verify status callbacks happened

		jobdirfile.recursivelyDeleteDirectory(); // clean up job dir
	}

	/**
	 * Test if a WADL document is available at the relative path
	 * "application.wadl".
	 */
	@Test
	public void testApplicationWadl() {
		WebResource webResource = resource();
		String serviceWadl = webResource.path("application.wadl")
				.accept(MediaTypes.WADL).get(String.class);

		assertTrue(serviceWadl.length() > 0);
	}
}
