package nl.nlesc.magma;

import java.io.IOException;
import java.net.URI;

import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.methods.HttpPut;
import org.apache.http.entity.StringEntity;
import org.apache.http.client.HttpClient;
import org.apache.http.impl.client.DefaultHttpClient;
import org.gridlab.gat.monitoring.MetricEvent;
import org.gridlab.gat.monitoring.MetricListener;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Listener on 'job.state' metric which writes state to file in job directory.
 *
 * @author Stefan Verhoeven <s.verhoeven@esciencecenter.nl>
 *
 */
public class JobStateListener implements MetricListener {
	protected final static Logger logger = LoggerFactory
			.getLogger(JobStateListener.class);
	protected URI status_cb_url;
	private HttpClient httpclient;

	/**
	 * Constructor
	 *
	 * @param status_cb_url
	 * @param httpclient
	 */
	public JobStateListener(URI status_cb_url) {
		this.httpclient = new DefaultHttpClient();
		this.status_cb_url = status_cb_url;
	}

	/**
	 * Constructor
	 *
	 * @param status_cb_url
	 * @param httpclient
	 */
	public JobStateListener(URI status_cb_url, HttpClient httpclient) {
		this.httpclient = httpclient;
		this.status_cb_url = status_cb_url;
	}

	/**
	 * Each time a event is fired the event.value is PUT to status_cb_url.
	 */
	@Override
	public void processMetricEvent(MetricEvent event) {
		try {
			HttpPut put = new HttpPut(status_cb_url);
			StringEntity body = new StringEntity(event.getValue().toString());
			put.setEntity(body);
			put.setHeader("Accept", "application/json");
			put.setHeader("Content-type", "text/plain; charset=UTF-8");
			httpclient.execute(put);
			put.releaseConnection();
		} catch (ClientProtocolException e) {
			logger.info("Unable to write job state to status callback url");
		} catch (IOException e) {
			logger.info("Unable to write job state to status callback url");
		}
		logger.debug("received state change: " + event.getValue() + " send to "
				+ status_cb_url);
	}

	public URI getStatusCallbackUrl() {
		return status_cb_url;
	}
}
