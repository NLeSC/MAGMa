package nl.esciencecenter.magma.mac;

import java.math.BigInteger;
import java.net.URI;
import java.security.InvalidKeyException;
import java.security.NoSuchAlgorithmException;
import java.security.SecureRandom;
import java.util.Date;
import java.util.Random;

import javax.crypto.Mac;
import javax.crypto.spec.SecretKeySpec;

import nl.esciencecenter.magma.gat.JobStateListener;

import org.apache.commons.codec.binary.Base64;
import org.apache.http.Header;
import org.apache.http.HttpRequest;
import org.apache.http.auth.AUTH;
import org.apache.http.auth.AuthenticationException;
import org.apache.http.auth.ContextAwareAuthScheme;
import org.apache.http.auth.Credentials;
import org.apache.http.auth.MalformedChallengeException;
import org.apache.http.client.methods.HttpUriRequest;
import org.apache.http.impl.client.RequestWrapper;
import org.apache.http.message.BasicHeader;
import org.apache.http.protocol.HttpContext;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class MacScheme implements ContextAwareAuthScheme {
	protected final static Logger logger = LoggerFactory
			.getLogger(JobStateListener.class);

	/** The name of this authorization scheme. */
	public static final String SCHEME_NAME = "MAC";
	private Date date = new Date();
	private Random random = new SecureRandom();

	public void processChallenge(Header header)
			throws MalformedChallengeException {
		// Not challenge based
	}

	public String getSchemeName() {
		return SCHEME_NAME;
	}

	public String getParameter(String name) {
		return null;
	}

	public String getRealm() {
		return null;
	}

	public boolean isConnectionBased() {
		return false;
	}

	public boolean isComplete() {
		return true;
	}

	@Deprecated
	public Header authenticate(Credentials credentials, HttpRequest request)
			throws AuthenticationException {
		return authenticate(credentials, request, null);
	}

	public Header authenticate(Credentials credentials, HttpRequest request,
			HttpContext context) throws AuthenticationException {
		String id = credentials.getUserPrincipal().getName();
		String key = credentials.getPassword();
		Long timestamp = getTimestamp();
		String nonce = getNonce();
		String data = getNormalizedRequestString((HttpUriRequest) request,
				nonce, timestamp);

		logger.warn(credentials.toString());

		String request_mac = calculateRFC2104HMAC(data, key,
				getAlgorithm(credentials));

		return new BasicHeader(AUTH.WWW_AUTH_RESP, headerValue(id, timestamp,
				nonce, request_mac));
	}

	private String headerValue(String id, Long timestamp, String nonce,
			String request_mac) {
		String headerValue = "MAC id=\"" + id + "\",";
		headerValue += "ts=\"" + timestamp + "\",";
		headerValue += "nonce=\"" + nonce + "\",";
		headerValue += "mac=\"" + request_mac + "\"";
		return headerValue;
	}

	/**
	 * Computes RFC 2104-compliant HMAC signature.
	 *
	 * @param data
	 *            The data to be signed.
	 * @param key
	 *            The signing key.
	 * @param algorithm
	 *            MAC algorithm implemented by javax.crypto.MAC
	 * @return The Base64-encoded RFC 2104-compliant HMAC signature.
	 * @throws RuntimeException
	 *             when signature generation fails
	 */
	private String calculateRFC2104HMAC(String data, String key,
			String algorithm) throws AuthenticationException {
		try {
			Mac mac = Mac.getInstance(algorithm);
			SecretKeySpec macKey = new SecretKeySpec(key.getBytes(), "RAW");
			mac.init(macKey);
			byte[] signature = mac.doFinal(data.getBytes());
			return Base64.encodeBase64String(signature);
		} catch (InvalidKeyException e) {
			throw new AuthenticationException("Failed to generate HMAC: "
					+ e.getMessage(), e);
		} catch (NoSuchAlgorithmException e) {
			throw new AuthenticationException("Algorithm is not supported", e);
		}
	}

	private String getNormalizedRequestString(HttpUriRequest request,
			String nonce, Long timestamp) {
		URI uri = request.getURI();
		// request can become wrapped, causing the request.getURI() to miss host
		// and port
		if (request instanceof RequestWrapper) {
			uri = ((HttpUriRequest) ((RequestWrapper) request).getOriginal())
					.getURI();
		}
		String normalized_request_string = timestamp + "\n";
		normalized_request_string += nonce + "\n";
		normalized_request_string += request.getMethod() + "\n";
		normalized_request_string += uri.getPath() + "\n";
		normalized_request_string += uri.getHost().toLowerCase() + "\n";
		normalized_request_string += getPort(uri) + "\n";
		normalized_request_string += "" + "\n";
		return normalized_request_string;
	}

	public Random getRandom() {
		return random;
	}

	public void setRandom(Random random) {
		this.random = random;
	}

	private String getNonce() {
		return new BigInteger(130, random).toString(32);
	}

	public void setDate(java.util.Date date) {
		this.date = date;
	}

	public Date getDate() {
		return date;
	}

	private Long getTimestamp() {
		Long timestamp = date.getTime() / 1000;
		return timestamp;
	}

	/**
	 *
	 * @param uri
	 * @return Port of `uri` based on explicit port or derived from scheme
	 */
	public static int getPort(URI uri) {
		int port = uri.getPort();
		if (port == -1) {
			String scheme = uri.getScheme();
			if (scheme.equals("http")) {
				port = 80;
			} else if (scheme.equals("https")) {
				port = 443;
			}
		}
		return port;
	}

	private String getAlgorithm(Credentials credentials) {
		String standardAlgo = ((MacCredential) credentials).getAlgorithm();
		return algorithmMapper(standardAlgo);
	}

	public static String algorithmMapper(String standard) {
		String java = standard;
		if (standard.equals("hmac-sha-1")) {
			java = "HmacSHA1";
		} else if (standard.equals("hmac-sha-256")) {
			java = "HmacSHA256";
		}
		// TODO implement registered extension algorithm
		return java;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		return true;
	}
}
