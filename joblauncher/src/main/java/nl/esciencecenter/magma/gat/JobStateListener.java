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

import com.google.common.base.Objects;

public class JobStateListener implements MetricListener {

	private URI status_callback_url;
	private HttpClient httpclient;

	public JobStateListener(URI status_callback_url, HttpClient httpclient) {
		this.status_callback_url = status_callback_url;
		this.httpclient = httpclient;
	}

	public void processMetricEvent(MetricEvent event) {
		String state = event.getValue().toString();
		HttpPut put = new HttpPut(status_callback_url);

		try {
			put.setEntity(new StringEntity(state));
			HttpResponse response = httpclient.execute(put);
		} catch (ClientProtocolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			put.releaseConnection();
		}
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		JobStateListener other = (JobStateListener) obj;
		return Objects.equal(this.status_callback_url, other.status_callback_url)
				&& Objects.equal(this.httpclient, other.httpclient);
	}

}
