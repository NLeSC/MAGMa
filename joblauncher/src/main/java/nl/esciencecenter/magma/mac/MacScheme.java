package nl.esciencecenter.magma.mac;

import java.net.URI;
import java.security.InvalidKeyException;
import java.security.NoSuchAlgorithmException;

import javax.crypto.Mac;
import javax.crypto.spec.SecretKeySpec;

import org.apache.commons.codec.binary.Base64;
import org.apache.http.Header;
import org.apache.http.HttpRequest;
import org.apache.http.auth.AuthScheme;
import org.apache.http.auth.AuthenticationException;
import org.apache.http.auth.Credentials;
import org.apache.http.auth.MalformedChallengeException;
import org.apache.http.client.methods.HttpUriRequest;
import org.apache.http.message.BasicHeader;

public class MacScheme implements AuthScheme {
	/** The name of this authorization scheme. */
	public static final String SCHEME_NAME = "MAC";

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

	public Header authenticate(Credentials credentials, HttpRequest request)
			throws AuthenticationException {
		String id = credentials.getUserPrincipal().getName();
		String key = credentials.getPassword();
		Long timestamp = getTimestamp();
		String nonce = getNonce();
		String data = getNormalizedRequestString((HttpUriRequest) request,
				nonce, timestamp);
		String request_mac = calculateRFC2104HMAC(data, key, getAlgorithm(credentials));

		return new BasicHeader("Authorization", headerValue(id, timestamp,
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
     * @param algorithm MAC algorithm implemented by javax.crypto.MAC
     * @return The Base64-encoded RFC 2104-compliant HMAC signature.
     * @throws RuntimeException
     *             when signature generation fails
     */
	private String calculateRFC2104HMAC(String data, String key, String algorithm) throws AuthenticationException {
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
			throw new AuthenticationException(algorithm
					+ " algorithm is not supported", e);
		}
	}

	private String getNormalizedRequestString(HttpUriRequest request,
			String nonce, Long timestamp) {
		URI uri = request.getURI();
		String normalized_request_string = timestamp + "\n";
		normalized_request_string += nonce + "\n";
		normalized_request_string += request.getMethod() + "\n";
		normalized_request_string += uri.getPath() + "\n";
		normalized_request_string += uri.getHost().toLowerCase() + "\n";
		normalized_request_string += getPort(uri) + "\n";
		normalized_request_string += "" + "\n";
		return normalized_request_string;
	}

	private String getNonce() {
		return java.util.UUID.randomUUID().toString();
	}

	private Long getTimestamp() {
		java.util.Date date = new java.util.Date();
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
		if (credentials instanceof MacCredential) {
			String standardAlgo = ((MacCredential) credentials).getAlgorithm();
			return algorithmMapper(standardAlgo);
		}
		return "HmacSHA1";
	}

	public static String algorithmMapper(String standard) {
		String java = "";
		if (standard == "hmac-sha-1") {
			java = "HmacSHA1";
		} else if (standard == "hmac-sha-256") {
			java = "HmacSHA256";
		}
		// TODO implement registered extension algorithm
		return java;
	}
}
