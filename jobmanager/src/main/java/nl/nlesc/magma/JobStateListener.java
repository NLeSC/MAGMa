package nl.nlesc.magma;

import java.io.IOException;
import java.io.InputStream;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.security.InvalidKeyException;
import java.security.NoSuchAlgorithmException;
import java.util.Properties;
import java.util.ResourceBundle;

import javax.crypto.Mac;
import javax.crypto.spec.SecretKeySpec;

import org.apache.commons.codec.binary.Base64;
import org.apache.http.HttpResponse;
import org.apache.http.client.ClientProtocolException;
import org.apache.http.client.HttpClient;
import org.apache.http.client.methods.HttpPut;
import org.apache.http.entity.StringEntity;
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
		String state = event.getValue().toString();

		putStatus(state);
	}

	public HttpResponse putStatus(String state) {
		// Credentials for MAGMa service
		// Retrieved from https://www.emetabolomics.org/magma/access_token
		// TODO put in config file

		String scheme = "MAC";
		String algorithm = "hmac-sha-1";
		String id = "some id";
		String key = "some key";
		String method = "PUT";

		String credentials = null;
		try {
			credentials = getCredentials(scheme, algorithm, id, key, method);
		} catch (InvalidKeyException e2) {
			logger.info(e2.toString());
			e2.printStackTrace();
		} catch (NoSuchAlgorithmException e2) {
			logger.info(e2.toString());
			e2.printStackTrace();
		}

		HttpPut put = new HttpPut(status_cb_url);
		put.setHeader("Accept", "application/json");
		put.setHeader("Content-type", "text/plain; charset=UTF-8");

		if (!credentials.isEmpty()) {
			put.setHeader("Authorization", credentials);
		}

		StringEntity body = null;
		try {
			body = new StringEntity(state);
		} catch (UnsupportedEncodingException e1) {
			logger.info(e1.toString());
			e1.printStackTrace();
		}
		put.setEntity(body);

		HttpResponse response = null;
		try {
			response = httpclient.execute(put);
		} catch (ClientProtocolException e) {
			logger.warn("Unable to write job state to status callback url");
			logger.debug(e.getMessage());
		} catch (IOException e) {
			logger.warn("Unable to write job state to status callback url");
			logger.debug(e.getMessage());
		} finally {
			put.releaseConnection();
		}
		logger.debug("Pushed state '" + state + "' to " + status_cb_url);
		return response;
	}

	/**
	 * @return Port of `status_cb_url` based on explicit port or derived from
	 *         scheme
	 */
	private int getPort() {
		int port = status_cb_url.getPort();
		if (port == -1) {
			String scheme = status_cb_url.getScheme();
			if (scheme.equals("http")) {
				port = 80;
			} else if (scheme.equals("https")) {
				port = 443;
			}
		}
		return port;
	}

	/**
	 * For now implements MAC Access Authentication See
	 * http://tools.ietf.org/html/draft-ietf-oauth-v2-http-mac-01
	 *
	 * @param scheme
	 *            "MAC" or returns empty string
	 * @param algorithm
	 *            the MAC algorithm
	 * @param id
	 * @param key
	 * @param method
	 * @return String for Authentication HTTP Header
	 * @throws NoSuchAlgorithmException
	 * @throws InvalidKeyException
	 */
	private String getCredentials(String scheme, String algorithm, String id,
			String key, String method) throws NoSuchAlgorithmException,
			InvalidKeyException {
		String credentials = "";

		if (scheme == "MAC") {
			String nonce = java.util.UUID.randomUUID().toString();
			java.util.Date date = new java.util.Date();
			Long timestamp = date.getTime() / 1000;

			String normalized_request_string = timestamp + "\n";
			normalized_request_string += nonce + "\n";
			normalized_request_string += method + "\n";
			normalized_request_string += status_cb_url.getPath() + "\n";
			normalized_request_string += status_cb_url.getHost().toLowerCase()
					+ "\n";
			normalized_request_string += getPort() + "\n";
			normalized_request_string += "" + "\n"; // ext

			// translate MAC algorithm from specification to Java Crypto
			// standard names
			if (algorithm == "hmac-sha-1") {
				algorithm = "HmacSHA1";
			} else if (algorithm == "hmac-sha-256") {
				algorithm = "HmacSHA256";
			}

			Mac mac = Mac.getInstance(algorithm);
			SecretKeySpec macKey = new SecretKeySpec(key.getBytes(), "RAW");
			mac.init(macKey);

			byte[] signature = mac
					.doFinal(normalized_request_string.getBytes());
			String request_mac = Base64.encodeBase64String(signature);
			credentials = "MAC id=\"" + id + "\",";
			credentials += "ts=\"" + timestamp + "\",";
			credentials += "nonce=\"" + nonce + "\",";
			credentials += "mac=\"" + request_mac + "\"";
		}
		return credentials;
	}

	public URI getStatusCallbackUrl() {
		return status_cb_url;
	}
}
