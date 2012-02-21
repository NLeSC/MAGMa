package org.nlesc.magma;

import java.io.FileWriter;
import java.io.IOException;

import org.gridlab.gat.io.File;
import org.gridlab.gat.monitoring.MetricEvent;
import org.gridlab.gat.monitoring.MetricListener;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Listener on 'job.state' metric which writes state to file in job directory.
 *
 * @author Stefan Verhoeven <s.verhoeven@esciencecenter.nl>
 *
 */
public class JobStateListener implements MetricListener {
    protected final static Logger logger = LoggerFactory
            .getLogger(JobStateListener.class);
    protected File jobdir;

    /**
     * Getter for job directory
     *
     * @return File
     */
    public File getJobdir() {
        return jobdir;
    }

    /**
     * Constructor
     *
     * @param jobdir
     * @throws IOException
     *             when jobdir is no directory
     */
    public JobStateListener(File jobdir) throws IOException {
        if (!jobdir.isDirectory()) {
            throw new IOException("Jobdir must be directory");
        }
        this.jobdir = jobdir;
    }

    /**
     * Each time a event is fired the event.value is written to jobdir/job.state file.
     */
    @Override
    public void processMetricEvent(MetricEvent event) {
        // write job status into file in jobdir so website can report it to user
        FileWriter fw;
        try {
            fw = new FileWriter(jobdir.getPath() + "/job.state");
            fw.write(event.getValue().toString());
            fw.close();
        } catch (IOException e) {
            logger.info("Unable to write job state to file");
        }
        logger.debug("received state change: " + event.getValue() + " of "
                + jobdir);
    }
}
