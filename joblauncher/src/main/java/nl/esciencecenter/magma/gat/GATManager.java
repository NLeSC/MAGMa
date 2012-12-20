package nl.esciencecenter.magma.gat;

import java.net.URISyntaxException;
import java.util.Set;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATContext;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.Preferences;
import org.gridlab.gat.URI;
import org.gridlab.gat.resources.ResourceBroker;

import com.yammer.dropwizard.lifecycle.Managed;

public class GATManager implements Managed {
	private final ResourceBroker broker;
	private final GATContext context;

	public GATManager(GATConfiguration configuration) throws GATObjectCreationException, URISyntaxException {
		context = GAT.getDefaultGATContext();

		// copy over preferences from config to default GAT context
    	Preferences prefs = new Preferences();
    	Set<String> keys = configuration.getPreferences().keySet();
    	for (String key : keys) {
            context.addPreference(key, configuration.getPreferences().get(key));
        }

		// TODO load SecurityContext from config, optionally interactive for passwords/passphrases

    	// create default broker
		URI brokerUri = new URI(configuration.getBrokerURI());
		broker = GAT.createResourceBroker(prefs, brokerUri);

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
