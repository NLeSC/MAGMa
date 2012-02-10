package org.nlesc.magma;


import java.net.URISyntaxException;

import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.resources.ResourceBroker;
import org.junit.Test;

import junit.framework.TestCase;

public class BrokerFactoryTest extends TestCase {

	@Test
	public void testGetBroker() throws GATObjectCreationException, URISyntaxException {
		BrokerFactory fact = new BrokerFactory();
		ResourceBroker broker = fact.getBroker();
		assertEquals(broker.toString(), "any://localhost");
	}

}
