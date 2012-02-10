/**
 *
 */
package org.nlesc.magma;

import java.net.URISyntaxException;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.URI;
import org.gridlab.gat.resources.ResourceBroker;

/**
 * @author stefanv
 * Toggles between local and glite broker
 *
 */
public class BrokerFactory {

	public ResourceBroker getBroker() throws GATObjectCreationException, URISyntaxException {
		return GAT.createResourceBroker(new URI(
                "any://localhost"));
	}
}
