package org.nlesc.magma;

import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.SoftwareDescription;
import org.nlesc.magma.entities.JobSubmitRequest;

public class JobDescriptionFactory {

	public JobDescription getJobDescription(JobSubmitRequest jobsubmission) throws Exception {
		if (jobsubmission.jobtype.equals("sleep")) {
				SoftwareDescription sd = new SoftwareDescription();
				// simulate job by sleeping a while
		        sd.setExecutable("/bin/sleep");
		        String[] args = { "30" };
		        sd.setArguments(args);

		        JobDescription jd = new JobDescription(sd);
		        return jd;
		} else {
			throw new Exception("Unknown job type: '" + jobsubmission.jobtype+"'");
		}
	}
}
