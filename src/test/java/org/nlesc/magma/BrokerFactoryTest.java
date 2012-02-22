package org.nlesc.magma;


import java.net.URISyntaxException;

import junit.framework.TestCase;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATContext;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.URI;
import org.gridlab.gat.resources.ResourceBroker;
import org.gridlab.gat.security.CertificateSecurityContext;
import org.junit.Test;

public class BrokerFactoryTest extends TestCase {

    @Test
    public void testGetBroker() throws GATObjectCreationException, URISyntaxException {
        BrokerFactory fact = new BrokerFactory();
        ResourceBroker broker = fact.getBroker();
        assertEquals(broker.toString(), "localq://localhost");
    }

    @Test
    public void testGetBrokerWithEmptyPassphrase() throws GATObjectCreationException, URISyntaxException {
        GATContext oldcontext = GAT.getDefaultGATContext();
        CertificateSecurityContext securityContext = new CertificateSecurityContext(
                new URI(System.getProperty("user.home") + "/.globus/userkey.pem"),
                new URI(System.getProperty("user.home") + "/.globus/usercert.pem"),
                "");

        GATContext context = new GATContext();
        context.addSecurityContext(securityContext);
        GAT.setDefaultGATContext(context);
        BrokerFactory fact = new BrokerFactory();

        ResourceBroker broker = fact.getBroker();
        assertEquals(broker.toString(), "localq://localhost");
        GAT.setDefaultGATContext(oldcontext);
    }

    @Test
    public void testGetBrokerWithPassphrase() throws GATObjectCreationException, URISyntaxException {
        GATContext oldcontext = GAT.getDefaultGATContext();
        CertificateSecurityContext securityContext = new CertificateSecurityContext(
                new URI(System.getProperty("user.home") + "/.globus/userkey.pem"),
                new URI(System.getProperty("user.home") + "/.globus/usercert.pem"),
                "blabla"); // a mocked password

        GATContext context = new GATContext();
        context.addSecurityContext(securityContext);
        GAT.setDefaultGATContext(context);
        BrokerFactory fact = new BrokerFactory();

        ResourceBroker broker = fact.getBroker();
        assertEquals(broker.toString(), "glite://wms4.grid.sara.nl:7443/glite_wms_wmproxy_server");
        GAT.setDefaultGATContext(oldcontext);
    }
}
