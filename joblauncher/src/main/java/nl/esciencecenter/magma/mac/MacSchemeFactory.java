package nl.esciencecenter.magma.mac;

import org.apache.http.auth.AuthScheme;
import org.apache.http.auth.AuthSchemeFactory;
import org.apache.http.params.HttpParams;

/**
 * Factory for {@Link MacScheme} implementations.
 *
 * @author verhoes
 */
public class MacSchemeFactory implements AuthSchemeFactory {
	/**
	 * Creates an instance of {@Link MacScheme}.
	 *
	 * @ return auth scheme
	 */
    public AuthScheme newInstance(HttpParams params) {
        return new MacScheme();
    }
}
