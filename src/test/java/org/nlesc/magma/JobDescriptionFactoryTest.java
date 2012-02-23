package org.nlesc.magma;

import static org.junit.Assert.assertArrayEquals;

import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

import junit.framework.TestCase;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.io.File;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.SoftwareDescription;
import org.junit.Test;
import org.nlesc.magma.entities.JobSubmitRequest;

public class JobDescriptionFactoryTest extends TestCase {
    protected JobDescriptionFactory fact;
    protected String jobdir;

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        jobdir = System.getProperty("java.io.tmpdir") + "/"
                + UUID.randomUUID().toString() + "/";
        fact = new JobDescriptionFactory();
    }

    @Override
    protected void tearDown() throws Exception {
        // clean up job dir
        GAT.createFile(jobdir).recursivelyDeleteDirectory();
        super.tearDown();
    }

    /**
     * Simplest job
     * @throws GATObjectCreationException
     */
    @Test
    public void testGetJobDescriptionHostname() throws GATObjectCreationException {
        JobSubmitRequest jobsubmission = new JobSubmitRequest();

        jobsubmission.jobdir = jobdir;
        jobsubmission.executable = "/bin/hostname";

        JobDescription jd = fact.getJobDescription(jobsubmission);
        SoftwareDescription sd = jd.getSoftwareDescription();

        assertEquals(sd.getExecutable(), "/bin/hostname");
        assertEquals(sd.getStdout().getAbsolutePath(), jobdir+"stdout.txt");
        assertEquals(sd.getStderr().getAbsolutePath(), jobdir+"stderr.txt");
        assertEquals(sd.getIntAttribute("time.max", 0), 0);
        assertEquals(sd.getIntAttribute("memory.min", 0), 0);
        assertEquals(sd.getIntAttribute("memory.max", 0), 0);
        String[] strarray = {};
        assertArrayEquals(sd.getArguments(), strarray);
        assertTrue(sd.getPostStaged().isEmpty());
        assertTrue(sd.getPreStaged().isEmpty());
    }

    /**
     * Job which uses all fields
     * @throws GATObjectCreationException
     */
    @Test
    public void testGetJobDescriptionWc() throws GATObjectCreationException {
        String[] args = { "-l", "etc/passwd", ">", "nr.accounts" };
        // use absolute file and relative file
        String[] prestage = { "/etc/passwd", "htaccess" };
        String[] poststage = { "nr.accounts" };

        JobSubmitRequest jobsubmission = new JobSubmitRequest(jobdir,
                "/bin/wc", args, prestage, poststage, "stderr.txt", "stdout.txt");

        jobsubmission.time_max = 60; // 1 hour
        jobsubmission.memory_min = 512; // at least 512Mb
        jobsubmission.memory_max = 2048; // at most 2048Mb

        JobDescription jd = fact.getJobDescription(jobsubmission);
        SoftwareDescription sd = jd.getSoftwareDescription();

        assertEquals(sd.getExecutable(), "/bin/wc");
        assertArrayEquals(sd.getArguments(), args);
        assertEquals(sd.getStderr().getAbsolutePath(), jobdir+"stderr.txt");
        assertEquals(sd.getStdout().getAbsolutePath(), jobdir+"stdout.txt");

        assertEquals(sd.getLongAttribute("time.max", 0), 60);
        assertEquals(sd.getIntAttribute("memory.min", 0), 512);
        assertEquals(sd.getIntAttribute("memory.max", 0), 2048);

        Map<File, File> actual_prestage = sd.getPreStaged();
        HashMap<File, File> expected_prestage = new HashMap<File, File>();
        expected_prestage.put(GAT.createFile("/etc/passwd"), null);
        expected_prestage.put(GAT.createFile(jobdir+"htaccess"), null);
        assertEquals(actual_prestage, expected_prestage);

        Map<File, File> actual_poststage = sd.getPostStaged();
        HashMap<File, File> expected_poststage = new HashMap<File, File>();
        expected_poststage.put(GAT.createFile("nr.accounts"), GAT.createFile(jobdir+"nr.accounts"));
        assertEquals(actual_poststage, expected_poststage);
    }
}
