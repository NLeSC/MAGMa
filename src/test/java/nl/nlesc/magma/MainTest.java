package nl.nlesc.magma;

import java.util.UUID;

import nl.nlesc.magma.Main;

import org.codehaus.jettison.json.JSONException;
import org.codehaus.jettison.json.JSONObject;
import org.glassfish.grizzly.http.server.HttpServer;
import org.gridlab.gat.GAT;
import org.gridlab.gat.GATInvocationException;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.io.File;

import com.sun.jersey.core.header.MediaTypes;
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
     * Test to see if a job can be submitted
     * This is a integration test, testing webservice and javagat local submitjob together
     * @throws JSONException
     * @throws GATObjectCreationException
     * @throws GATInvocationException
     */
    public void testJobResourcePOST() throws JSONException, GATObjectCreationException, GATInvocationException {
        // states of job are written to jobdir, so need an temp dir
        String jobdir = System.getProperty("java.io.tmpdir")+"/"+UUID.randomUUID().toString();
        File jobdirfile = GAT.createFile(jobdir);
        jobdirfile.mkdir();

        JSONObject requestMsg = new JSONObject();
        String[] args = { "1" };
        requestMsg.put("jobdir", jobdir).put("executable", "sleep").put("arguments", args);

        JSONObject responseMsg = r.path("job").post(JSONObject.class, requestMsg);

        assertEquals("0", responseMsg.get("jobid")); // first job == 0 and GAT is restarted each time

        GAT.end(); // stop job
        jobdirfile.recursivelyDeleteDirectory(); // clean up job dir with state file written by state listener
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
