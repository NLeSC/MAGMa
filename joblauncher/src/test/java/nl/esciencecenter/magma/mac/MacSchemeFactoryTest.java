package nl.esciencecenter.magma.mac;

import static org.junit.Assert.*;

import org.apache.http.auth.AuthScheme;
import org.apache.http.params.BasicHttpParams;
import org.apache.http.params.HttpParams;
import org.junit.Test;

public class MacSchemeFactoryTest {

	@Test
	public void testNewInstance() {
		MacSchemeFactory factory = new MacSchemeFactory();
		HttpParams params = new BasicHttpParams();
		AuthScheme scheme = factory.newInstance(params);

		assertEquals(new MacScheme(), scheme);
	}

}
