package nl.nlesc.magma;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.util.UUID;

import junit.framework.TestCase;

import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpPut;
import org.gridlab.gat.monitoring.MetricEvent;
import org.gridlab.gat.resources.Job.JobState;
import org.junit.Test;
import org.mockito.ArgumentCaptor;

public class JobStateListenerTest extends TestCase {
	private String jobid;
	private URI status_cb_url;

	@Override
	protected void setUp() throws Exception {
		super.setUp();
		jobid = UUID.randomUUID().toString();
		status_cb_url = new URI("http://example.com/status/" + jobid + ".json");
	}

	@Override
	protected void tearDown() throws Exception {
		super.tearDown();
	}

	@Test
	public void testJobStateListener() throws IOException {
		JobStateListener listener = new JobStateListener(status_cb_url);
		assertEquals(status_cb_url, listener.getStatusCallbackUrl());
	}

	@Test
	public void testProcessMetricEvent() throws IOException {
		HttpClient ua = mock(HttpClient.class);
		JobStateListener listener = new JobStateListener(status_cb_url, ua);
		// cannot equal the request, so have to capture
		ArgumentCaptor<HttpPut> httpRequest = ArgumentCaptor
				.forClass(HttpPut.class);

		MetricEvent event = mock(MetricEvent.class);
		when(event.getValue()).thenReturn(JobState.STOPPED);

		listener.processMetricEvent(event);

		verify(ua).execute(httpRequest.capture());
		assertEquals(httpRequest.getValue().getMethod(), "PUT");
		assertEquals(httpRequest.getValue().getURI(), status_cb_url);
		assertEquals(getRequestBody(httpRequest), JobState.STOPPED.toString());
	}

	private String getRequestBody(ArgumentCaptor<HttpPut> httpRequest)
			throws IOException {
		InputStream input = httpRequest.getValue().getEntity().getContent();
		StringBuilder builder = new StringBuilder();

		int data = input.read();
		while (data != -1) {
			builder.append((char) data);
			data = input.read();
		}
		String body = builder.toString();
		return body;
	}
}
