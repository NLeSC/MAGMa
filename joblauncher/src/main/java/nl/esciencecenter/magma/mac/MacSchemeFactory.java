package nl.esciencecenter.magma.mac;

import org.apache.http.auth.AuthScheme;
import org.apache.http.auth.AuthSchemeFactory;
import org.apache.http.params.HttpParams;

public class MacSchemeFactory implements AuthSchemeFactory {
	public AuthScheme newInstance(HttpParams params) {
		return new MacScheme();
	}
}
