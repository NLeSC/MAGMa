package nl.esciencecenter.magma.mac;

import static com.yammer.dropwizard.testing.JsonHelpers.fromJson;
import static com.yammer.dropwizard.testing.JsonHelpers.jsonFixture;
import static org.fest.assertions.api.Assertions.assertThat;
import static org.hamcrest.Matchers.is;
import static org.junit.Assert.assertThat;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;

import org.apache.http.auth.AuthScope;
import org.apache.http.auth.BasicUserPrincipal;
import org.junit.Before;
import org.junit.Test;

public class MacCredentialTest {
	MacCredential cred;

	@Before
	public void setup() {
		try {
			cred = new MacCredential("id", "key", new URI("http://localhost"));
		} catch (URISyntaxException e) {
			e.printStackTrace();
		}
	}

	@Test
	public void testMacCredentialStringStringURI() throws URISyntaxException {
		assertThat(cred.getAlgorithm()).isEqualTo("hmac-sha-1");
		assertThat(cred.getId()).isEqualTo("id");
		assertThat(cred.getKey()).isEqualTo("key");
		assertThat(cred.getScope()).isEqualTo(new URI("http://localhost"));
	}

	@Test
	public void testSetId() {
		cred.setId("id2");
		assertThat(cred.getId()).isEqualTo("id2");
	}

	@Test
	public void testSetKey() {
		cred.setKey("key2");
		assertThat(cred.getKey()).isEqualTo("key2");
	}

	@Test
	public void testSetScope() throws URISyntaxException {
		URI uri = new URI("https://example.com");
		cred.setScope(uri);
		assertThat(cred.getScope()).isEqualTo(uri);
	}

	@Test
	public void testSetAlgorithm() {
		cred.setAlgorithm("hmac-sha-256");
		assertThat(cred.getAlgorithm()).isEqualTo("hmac-sha-256");
	}

	@Test
	public void testGetAuthScope() {
		AuthScope expected = new AuthScope("localhost", 80, "", "MAC");
		assertThat(cred.getAuthScope()).isEqualTo(expected);
	}

	@Test
	public void testGetUserPrincipal() {
		assertThat(cred.getUserPrincipal()).isEqualTo(
				new BasicUserPrincipal("id"));
	}

	@Test
	public void testGetPassword() {
		assertThat(cred.getPassword()).isEqualTo("key");

	}

	@Test
	public void deserializesFromJSON() throws IOException {
		assertThat(
				"a MacCredential can be deserialized from JSON",
				fromJson(jsonFixture("fixtures/mac_credential.json"),
						MacCredential.class), is(cred));
	}

	@Test
	public void hasAWorkingEqualsMethod() throws Exception {
		assertThat(cred.equals(cred)).isTrue();

		assertThat(
				cred.equals(new MacCredential("id", "key", new URI(
						"http://localhost")))).isTrue();

		assertThat(cred.equals(null)).isFalse();

		assertThat(cred.equals("string")).isFalse();

		assertThat(
				cred.equals(new MacCredential("id2", "key", new URI(
						"http://localhost")))).isFalse();

		assertThat(
				cred.equals(new MacCredential("id", "key2", new URI(
						"http://localhost")))).isFalse();

		assertThat(
				cred.equals(new MacCredential("id2", "key", new URI(
						"https://example.com")))).isFalse();
	}

	@Test
	public void testToString() {
		assertThat(cred.toString()).isEqualTo("MacCredential [id=id, key=key, algorithm=hmac-sha-1, scope=http://localhost]");
	}
}
