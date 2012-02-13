package org.nlesc.magma;


import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.when;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.UUID;

import junit.framework.TestCase;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.io.File;
import org.gridlab.gat.monitoring.MetricEvent;
import org.gridlab.gat.resources.Job.JobState;
import org.junit.Test;

public class JobStateListenerTest extends TestCase {
	private String jobdir;
	private File jobdirfile;

	@Override
	protected void setUp() throws Exception {
		super.setUp();
		jobdir = System.getProperty("java.io.tmpdir")+"/"+UUID.randomUUID().toString();
		jobdirfile = GAT.createFile(jobdir);
		jobdirfile.mkdir();
	}

	@Override
	protected void tearDown() throws Exception {
		jobdirfile.recursivelyDeleteDirectory();
		super.tearDown();
	}

	@Test
	public void testJobStateListenerNoDir() throws GATObjectCreationException {
		File f = GAT.createFile("bla");
		try {
			new JobStateListener(f);
			fail("Job state listener requires a directory");
		} catch (IOException e) {
			assertTrue(true);
		} finally {
			f.delete();
		}
	}

	@Test
	public void testJobStateListener() throws IOException {
		JobStateListener listener = new JobStateListener(jobdirfile);
		assertEquals(jobdirfile, listener.getJobdir());
	}

	@Test
	public void testProcessMetricEventStopped() throws IOException {
		JobStateListener listener = new JobStateListener(jobdirfile);

		MetricEvent event = mock(MetricEvent.class);
		when(event.getValue()).thenReturn(JobState.STOPPED);

		listener.processMetricEvent(event);

		FileReader fr = new FileReader(jobdir+"/job.state");
		BufferedReader br = new BufferedReader(fr);
		String state = br.readLine();
		fr.close();
		assertEquals(JobState.STOPPED.toString(), state);
	}
}
