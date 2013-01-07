package nl.esciencecenter.magma.gat;

import java.io.IOException;
import java.net.URI;

import org.apache.http.HttpResponse;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpPut;
import org.apache.http.entity.StringEntity;
import org.gridlab.gat.monitoring.MetricEvent;
import org.gridlab.gat.monitoring.MetricListener;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.base.Objects;

/**
 * Listener for job state
 *
 * @author verhoes
 *
 */
public class JobStateListener implements MetricListener {
    protected final static Logger logger = LoggerFactory
            .getLogger(JobStateListener.class);
    private URI status_callback_url;
    private HttpClient httpclient;

    /**
     *
     * @param status_callback_url URL to put job state to
     * @param httpclient Client to put job state to status_callback_url
     */
    public JobStateListener(URI status_callback_url, HttpClient httpclient) {
        this.status_callback_url = status_callback_url;
        this.httpclient = httpclient;
    }

    /**
     * PUTs event value to status_callback_url using httpclient
     */
    public void processMetricEvent(MetricEvent event) {
        String state = event.getValue().toString();
        HttpPut put = new HttpPut(status_callback_url);

        try {
            put.setEntity(new StringEntity(state));
            HttpResponse response = httpclient.execute(put);
            logger.info("Send '" + state + "' to " + status_callback_url
                    + " returned " + response.getStatusLine());
        } catch (ClientProtocolException e) {
            logger.warn("Unable to write job state to status callback url", e);
        } catch (IOException e) {
            logger.warn("Unable to write job state to status callback url", e);
        } finally {
            put.releaseConnection();
        }
    }

    @Override
    public int hashCode() {
        return Objects.hashCode(status_callback_url, httpclient);
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        JobStateListener other = (JobStateListener) obj;
        return Objects.equal(this.status_callback_url,
                other.status_callback_url)
                && Objects.equal(this.httpclient, other.httpclient);
    }

}
