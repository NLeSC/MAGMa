package nl.esciencecenter.magma.api;

import static org.junit.Assert.*;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.HashMap;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.io.File;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.SoftwareDescription;
import org.junit.Before;
import org.junit.Test;

import static com.yammer.dropwizard.testing.JsonHelpers.*;
import static org.fest.assertions.api.Assertions.assertThat;
import static org.hamcrest.Matchers.*;

public class JobSubmitRequestTest {
    private JobSubmitRequest request;

    @Before
    public void setUp() {
        request = sampleRequest();
    }

    @Test
    public void testJobSubmitRequestStringStringStringArrayStringArrayStringArrayStringStringURI() {
        assertNotNull(request);
    }

    @Test
    public void testGetStatus_callback_url() throws URISyntaxException {
        assertEquals(request.getStatus_callback_url(), new URI(
                "http://localhost/status"));
    }

    @Test
    public void testJobSubmitRequest() {
        assertNotNull(new JobSubmitRequest());
    }

    @Test
    public void testToJobDescription() throws GATObjectCreationException {
        JobDescription jd = request.toJobDescription();

        SoftwareDescription sd = new SoftwareDescription();
        sd.setExecutable("/bin/sh");
        String[] arguments = { "magma.sh" };
        sd.setArguments(arguments);
        sd.setStderr(GAT.createFile("/tmp/jobdir/stderr.txt"));
        sd.setStdout(GAT.createFile("/tmp/jobdir/stdout.txt"));
        sd.addPreStagedFile(GAT.createFile("/tmp/jobdir/magma.sh"));
        sd.addPreStagedFile(GAT.createFile("/tmp/jobdir/data.mzxml"));
        sd.addPostStagedFile(GAT.createFile("results.db"), GAT.createFile("/tmp/jobdir/results.db"));
        JobDescription ejd =  new JobDescription(sd);
        assertEquals(jd.toString(), ejd.toString());
    }

    @Test
    public void testToJobDescriptionAttributes() throws GATObjectCreationException {
        request.time_max = 60; // 1 hour
        request.memory_min = 512; // at least 512Mb
        request.memory_max = 2048; // at most 2048Mb

        JobDescription jd = request.toJobDescription();

        SoftwareDescription sd = jd.getSoftwareDescription();
        assertEquals(sd.getLongAttribute("time.max", 0), 60);
        assertEquals(sd.getIntAttribute("memory.min", 0), 512);
        assertEquals(sd.getIntAttribute("memory.max", 0), 2048);
    }

    @Test
    public void testToJobDescriptionAbsolutePreStage() throws GATObjectCreationException {
        String[] prestaged = {"/etc/passwd"};
        request.prestaged = prestaged;

        JobDescription jd = request.toJobDescription();

        HashMap<File, File> expected_prestage = new HashMap<File, File>();
        expected_prestage.put(GAT.createFile("/etc/passwd"), null);
        SoftwareDescription sd = jd.getSoftwareDescription();
        assertEquals(expected_prestage, sd.getPreStaged());
    }

    private JobSubmitRequest sampleRequest() {
        String[] arguments = { "magma.sh" };
        String[] prestaged = { "magma.sh", "data.mzxml" };
        String[] poststaged = { "results.db" };
        URI cb = null;
        try {
            cb = new URI("http://localhost/status");
        } catch (URISyntaxException e) {
        }
        return new JobSubmitRequest("/tmp/jobdir/", "/bin/sh", arguments,
                prestaged, poststaged, "stderr.txt", "stdout.txt", cb);
    }

    @Test
    public void serializesToJSON() throws IOException {
        assertThat("a JobSubmitRequest can be serialized to JSON",
                   asJson(request),
                   is(equalTo(jsonFixture("fixtures/request.json"))));
    }

    @Test
    public void deserializesFromJSON() throws IOException {
        assertThat(
                "a JobSubmitRequest can be deserialized from JSON",
                fromJson(jsonFixture("fixtures/request.json"),
                        JobSubmitRequest.class), is(request));
    }

    @Test
    public void testEquals() {
        assertThat(request.equals(null)).isFalse();

        assertThat(request.equals("string")).isFalse();

        assertTrue(request.equals(request));

        assertTrue(request.equals(sampleRequest()));

        JobSubmitRequest r2 = sampleRequest();
        r2.executable = "/bin/bash";
        assertFalse(request.equals(r2));

        JobSubmitRequest r3 = sampleRequest();
        r3.jobdir = "/tmp/jobdir2";
        assertFalse(request.equals(r3));

        JobSubmitRequest r4 = sampleRequest();
        r4.stderr = "error.log";
        assertFalse(request.equals(r4));

        JobSubmitRequest r5 = sampleRequest();
        r5.stdout = "out.log";
        assertFalse(request.equals(r5));

        JobSubmitRequest r6 = sampleRequest();
        r6.memory_max = 512;
        assertFalse(request.equals(r6));

        JobSubmitRequest r7 = sampleRequest();
        r7.memory_min = 512;
        assertFalse(request.equals(r7));

        JobSubmitRequest r8 = sampleRequest();
        r8.time_max = 60;
        assertFalse(request.equals(r8));

        try {
            JobSubmitRequest r9 = sampleRequest();
            r9.status_callback_url = new URI("http://example.com");
            assertFalse(request.equals(r9));
        } catch (URISyntaxException e) {
            fail();
        }
    }

    @Test
    public void testToString() {
        String s = "JobSubmitRequest{jobdir=/tmp/jobdir/, executable=/bin/sh, stderr=stderr.txt, stdout=stdout.txt, arguments=[magma.sh], prestaged=[magma.sh, data.mzxml], poststaged=[results.db], time_max=0, memory_min=0, memory_max=0, status_callback_url=http://localhost/status}";
        assertEquals(s, request.toString());
    }
}