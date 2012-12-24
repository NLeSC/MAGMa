package nl.esciencecenter.magma.gat;

import static org.fest.assertions.api.Assertions.assertThat;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;
import static org.mockito.Matchers.any;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;
import static org.mockito.Mockito.when;

import java.io.IOException;
import java.io.InputStreamReader;
import java.net.URI;
import java.net.URISyntaxException;

import org.apache.http.HttpResponse;
import org.apache.http.HttpVersion;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpPut;
import org.apache.http.client.methods.HttpUriRequest;
import org.apache.http.impl.client.DefaultHttpClient;
import org.apache.http.message.BasicStatusLine;
import org.gridlab.gat.monitoring.MetricEvent;
import org.gridlab.gat.resources.Job.JobState;
import org.junit.Test;
import org.mockito.ArgumentCaptor;

import com.google.common.io.CharStreams;

public class JobStateListenerTest {

    @Test
    public void testJobStateListener() throws URISyntaxException {
        URI cb = new URI("http://localhost/status");
        HttpClient ua = new DefaultHttpClient();
        assertNotNull(new JobStateListener(cb, ua));
    }

    @Test
    public void testProcessMetricEvent() throws URISyntaxException,
            ClientProtocolException, IOException {
        URI cb = new URI("http://localhost/status");
        HttpClient ua = mock(HttpClient.class);
        HttpResponse response = mock(HttpResponse.class);
        when(response.getStatusLine()).thenReturn(
                new BasicStatusLine(HttpVersion.HTTP_1_1, 200, "OK"));
        when(ua.execute(any(HttpUriRequest.class))).thenReturn(response);
        JobStateListener l = new JobStateListener(cb, ua);
        MetricEvent event = mock(MetricEvent.class);
        when(event.getValue()).thenReturn(JobState.STOPPED);

        l.processMetricEvent(event);

        // HttpPut.equals() fails when creating put object here, so capture it
        // to verify
        ArgumentCaptor<HttpPut> argument = ArgumentCaptor
                .forClass(HttpPut.class);
        verify(ua).execute(argument.capture());
        HttpPut r = argument.getValue();
        assertEquals("PUT", r.getMethod());
        assertEquals("http://localhost/status", r.getURI().toString());
        String content = CharStreams.toString(new InputStreamReader(r
                .getEntity().getContent()));
        assertEquals(JobState.STOPPED.toString(), content);
    }

    @Test
    public void testEquals() throws URISyntaxException {
        URI cb = new URI("http://localhost/status");
        HttpClient ua = new DefaultHttpClient();
        JobStateListener l = new JobStateListener(cb, ua);

        assertTrue(l.equals(l));

        assertThat(l.equals(null)).isFalse();

        assertThat(l.equals("string")).isFalse();

        URI cb2 = new URI("https://example.com/status");
        JobStateListener l2 = new JobStateListener(cb2, ua);
        assertFalse(l.equals(l2));

        HttpClient ua2 = new DefaultHttpClient();
        ua2.getParams().setParameter("key", "val");
        JobStateListener l3 = new JobStateListener(cb, ua2);
        assertFalse(l.equals(l3));
    }
}
