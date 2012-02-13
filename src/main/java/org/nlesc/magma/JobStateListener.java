package org.nlesc.magma;

import java.io.FileWriter;
import java.io.IOException;

import org.gridlab.gat.io.File;
import org.gridlab.gat.monitoring.MetricEvent;
import org.gridlab.gat.monitoring.MetricListener;

public class JobStateListener implements MetricListener {
	private File jobdir;

	public File getJobdir() {
		return jobdir;
	}

	public JobStateListener(File jobdir) throws IOException {
		if (!jobdir.isDirectory()) {
			throw new IOException("Jobdir must be directory");
		}
		this.jobdir = jobdir;
	}

	public void processMetricEvent(MetricEvent event) {
		// write job status into file in jobdir so website can report it to user
		FileWriter fw;
		try {
			fw = new FileWriter(jobdir.getPath()+"/job.state");
			fw.write(event.getValue().toString());
			fw.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		System.err.println("received state change: " + event.getValue() + " of " + jobdir);
	}
}
