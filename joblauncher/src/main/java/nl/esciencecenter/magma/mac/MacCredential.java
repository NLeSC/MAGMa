package nl.esciencecenter.magma.mac;

import java.net.URI;
import java.security.Principal;

import javax.validation.constraints.NotNull;

import org.apache.http.auth.AuthScope;
import org.apache.http.auth.BasicUserPrincipal;
import org.apache.http.auth.Credentials;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.google.common.base.Objects;

public class MacCredential implements Credentials {
	@NotNull
	@JsonProperty
	private String id = null;

	@NotNull
	@JsonProperty
	private String key = null;

	@JsonProperty
	@NotNull
	private String algorithm = "hmac-sha-1";

	@JsonProperty
	@NotNull
	private URI scope = null;

	public MacCredential(String id, String key, URI scope) {
		this.id = id;
		this.key = key;
		this.scope = scope;
	}

	public MacCredential() {
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getKey() {
		return key;
	}

	public void setKey(String key) {
		this.key = key;
	}

	public URI getScope() {
		return scope;
	}

	public void setScope(URI scope) {
		this.scope = scope;
	}

	public String getAlgorithm() {
		return algorithm;
	}

	public void setAlgorithm(String algorithm) {
		this.algorithm = algorithm;
	}

	public AuthScope getAuthScope() {
		return new AuthScope(scope.getHost(), MacScheme.getPort(scope), "",
				MacScheme.SCHEME_NAME);
	}

	public Principal getUserPrincipal() {
		return new BasicUserPrincipal(id);
	}

	public String getPassword() {
		return key;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		MacCredential other = (MacCredential) obj;
		return Objects.equal(this.id, other.id)
				&& Objects.equal(this.key, other.key)
				&& Objects.equal(this.algorithm, other.algorithm)
				&& Objects.equal(this.scope, other.scope);
	}
}
