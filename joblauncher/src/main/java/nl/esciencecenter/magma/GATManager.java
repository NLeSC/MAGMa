package nl.esciencecenter.magma;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATContext;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.resources.ResourceBroker;

import com.yammer.dropwizard.lifecycle.Managed;

public class GATManager implements Managed {
	private final ResourceBroker broker;
	private final GATContext context;

	public GATManager(GATConfiguration configuration) throws GATObjectCreationException {
		broker = GAT.createResourceBroker(configuration.getBrokerURI());
		context = new GATContext();
		context.addPreferences(configuration.getPreferences());
		GAT.setDefaultGATContext(context);
	}

	public void start() throws Exception {
	}

	public void stop() throws Exception {
		GAT.end();
	}

	public ResourceBroker getBroker() {
		return broker;
	}
}
