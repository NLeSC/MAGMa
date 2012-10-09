package nl.nlesc.magma.resources;

import static org.mockito.Matchers.any;
import static org.mockito.Matchers.anyObject;
import static org.mockito.Matchers.anyString;
import static org.mockito.Matchers.eq;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import java.util.UUID;

import junit.framework.TestCase;

import nl.nlesc.magma.BrokerFactory;
import nl.nlesc.magma.JobStateListener;
import nl.nlesc.magma.entities.JobSubmitRequest;
import nl.nlesc.magma.entities.JobSubmitResponse;
import nl.nlesc.magma.resources.JobResource;

import org.gridlab.gat.GAT;
import org.gridlab.gat.io.File;
import org.gridlab.gat.monitoring.MetricListener;
import org.gridlab.gat.resources.Job;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.ResourceBroker;
import org.junit.Test;

public class JobResourceTest extends TestCase {
    JobResource resource;
    ResourceBroker mockedbroker;

    @Override
    protected void setUp() throws Exception {
        super.setUp();
        resource = new JobResource();
        // mock getBroker() so no job is actually run
        BrokerFactory mockedbrokerfactory = mock(BrokerFactory.class);
        mockedbroker = mock(ResourceBroker.class);
        when(mockedbrokerfactory.getBroker()).thenReturn(mockedbroker);
        resource.setBroker(mockedbrokerfactory);
    }

    @Override
    protected void tearDown() throws Exception {
        GAT.end();
        super.tearDown();
    }

    @Test
    public void testSubmitJob() throws Exception {
        Job mockedjob = mock(Job.class);
        when(mockedjob.getJobID()).thenReturn(12345);
        when(
                mockedbroker.submitJob(any(JobDescription.class),
                        (MetricListener) anyObject(), anyString())).thenReturn(
                mockedjob);
        File jobdirfile = GAT.createFile(System.getProperty("java.io.tmpdir")
                + "/" + UUID.randomUUID().toString());
        jobdirfile.mkdir();
        String[] args = { "1" };
        String[] stage = {};
        JobSubmitRequest in = new JobSubmitRequest(jobdirfile.getPath(),
                "sleep", args, stage, stage, "stderr.txt", "stdout.txt");

        JobSubmitResponse out = resource.submitJob(in);

        assertEquals("12345", out.jobid);
        verify(mockedbroker).submitJob(any(JobDescription.class),
                any(JobStateListener.class), eq("job.status"));
        GAT.end();
        jobdirfile.recursivelyDeleteDirectory();
    }
}
