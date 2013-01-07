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
    public AuthScheme newInstance(final HttpParams params) {
        return new MacScheme();
    }
}
