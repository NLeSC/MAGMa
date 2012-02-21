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

	@Test
	public void testGetJobDescriptionWc() throws GATObjectCreationException {
		String[] args = { "-l", "etc/passwd", ">", "nr.accounts" };
		// use absolute file and relative file
		String[] prestage = { "/etc/passwd", "htaccess" };
		String[] poststage = { "nr.accounts" };

		JobSubmitRequest jobsubmission = new JobSubmitRequest(jobdir,
                "/bin/wc", args, prestage, poststage, "stderr.txt", "stdout.txt");

		JobDescription jd = fact.getJobDescription(jobsubmission);
		SoftwareDescription sd = jd.getSoftwareDescription();

        assertEquals(sd.getExecutable(), "/bin/wc");
        assertArrayEquals(sd.getArguments(), args);
        assertEquals(sd.getStderr().getAbsolutePath(), jobdir+"stderr.txt");
        assertEquals(sd.getStdout().getAbsolutePath(), jobdir+"stdout.txt");

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
