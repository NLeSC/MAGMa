package org.nlesc.magma;

import org.gridlab.gat.monitoring.MetricEvent;
import org.gridlab.gat.monitoring.MetricListener;
import org.gridlab.gat.resources.Job.JobState;

public class JobSubmitCallback implements MetricListener {
	private String jobdir;

	public JobSubmitCallback(String jobdir) {
		this.jobdir = jobdir;
	}

	public void processMetricEvent(MetricEvent event) {
		System.err.println("received state change: " + event.getValue() + " of " + jobdir);
		// TODO write job status into file in jobdir so website can report it to user
        if (event.getValue() == JobState.STOPPED) {
    		System.err.println("Job completed");
    		// TODO copy results files back to jobdir
        }
	}
}
