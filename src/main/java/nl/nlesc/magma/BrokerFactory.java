/**
 *
 */
package nl.nlesc.magma;

import java.net.URISyntaxException;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.URI;
import org.gridlab.gat.resources.ResourceBroker;

/**
 * Factory to create resource brokers.
 *
 * @author Stefan Verhoeven <s.verhoeven@esciencecenter.nl>
 *
 */
public class BrokerFactory {

    /**
     * Creates resource broker. If passphrase was given during startup returns
     * glite broker on grid.sara.nl else returns localhost broker.
     *
     * @return ResourceBroker
     * @throws GATObjectCreationException
     * @throws URISyntaxException
     */
    public ResourceBroker getBroker() throws GATObjectCreationException,
            URISyntaxException {
        String brokeruri = "localq://localhost";
        if (!GAT.getDefaultGATContext().getSecurityContexts().isEmpty()
                && GAT.getDefaultGATContext().getSecurityContexts().get(0)
                        .getPassword().length() > 0) {
            brokeruri = "glite://wms4.grid.sara.nl:7443/glite_wms_wmproxy_server";
        }
        return GAT.createResourceBroker(new URI(brokeruri));
    }
}
