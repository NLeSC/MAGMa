package nl.esciencecenter.magma;

import java.util.Set;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATContext;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.Preferences;
import org.gridlab.gat.resources.ResourceBroker;

import com.yammer.dropwizard.lifecycle.Managed;

public class GATManager implements Managed {
	private final ResourceBroker broker;
	private final GATContext context;

	public GATManager(GATConfiguration configuration) throws GATObjectCreationException {
		context = new GATContext();

		// copy over preferences
    	Preferences prefs = new Preferences();
    	Set<String> keys = configuration.getPreferences().keySet();
    	for (String key : keys) {
            prefs.put(key, configuration.getPreferences().get(key));
        }
		context.addPreferences(prefs);

		broker = GAT.createResourceBroker(context, configuration.getBrokerURI());
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
