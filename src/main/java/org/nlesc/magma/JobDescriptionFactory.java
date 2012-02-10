package org.nlesc.magma;

import org.gridlab.gat.GAT;
import org.gridlab.gat.Preferences;
import org.gridlab.gat.URI;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.SoftwareDescription;
import org.gridlab.gat.resources.WrapperJobDescription;
import org.gridlab.gat.resources.WrapperSoftwareDescription;
import org.nlesc.magma.entities.JobSubmitRequest;

public class JobDescriptionFactory {

	/**
	 *
	 * @param jobsubmission
	 * @return
	 * @throws Exception
	 */
	public JobDescription getJobDescription(JobSubmitRequest jobsubmission) throws Exception {
		if (jobsubmission.jobtype.equals("sleep")) {
				SoftwareDescription sd = new SoftwareDescription();
				// simulate job by sleeping a while
		        sd.setExecutable("/bin/sleep");
		        String[] args = { "30" };
		        sd.setArguments(args);

		        JobDescription jd = new JobDescription(sd);
		        return jd;
		} else if (jobsubmission.jobtype.equals("mzxmllocal")) {

			SoftwareDescription sd = new SoftwareDescription();
			// simulate job by sleeping a while
	        sd.setExecutable("/usr/bin/env");
	        String [] args = {
	        		"mscore_mzxml", // expects mscore_mzxml to be in path
	        		// "echo",
					"-p", jobsubmission.arguments.precision,
					"-c", jobsubmission.arguments.mscutoff,
					"-d", jobsubmission.arguments.msmscutoff,
					"-i", jobsubmission.arguments.ionisation,
					"-n", jobsubmission.arguments.nsteps,
					"-m", jobsubmission.arguments.phase,
					// input files
					"-f", jobsubmission.jobdir+"/data.mzxml",
					"-s", jobsubmission.jobdir+"/smiles.sd",
					// output file
					"-o", jobsubmission.jobdir+"/results.db"
			};
	        sd.setArguments(args);

	        sd.setStdout(GAT.createFile(jobsubmission.jobdir+"/stdout.txt"));
	        sd.setStderr(GAT.createFile(jobsubmission.jobdir+"/stderr.txt"));

	        JobDescription jd = new JobDescription(sd);
	        return jd;
		} else if (jobsubmission.jobtype.equals("mzxmlremote")) {
	        WrapperSoftwareDescription wsd = new WrapperSoftwareDescription();
	        wsd.setStdout(GAT.createFile("wrapper.stdout"));
	        wsd.setStderr(GAT.createFile("wrapper.stderr"));
	        wsd.setExecutable("/usr/local/package/jdk1.6.0/bin/java");

	        WrapperJobDescription wjd = new WrapperJobDescription(wsd);

	        SoftwareDescription sd = new SoftwareDescription();
	        // TODO stage magma program tarball
	        // TODO replace args[0] with staged shell script
	        sd.setExecutable("/usr/bin/env");
	        String [] args = {
	        		"mscore_mzxml", // expects mscore_mzxml to be in path
	        		// "echo",
					"-p", jobsubmission.arguments.precision,
					"-c", jobsubmission.arguments.mscutoff,
					"-d", jobsubmission.arguments.msmscutoff,
					"-i", jobsubmission.arguments.ionisation,
					"-n", jobsubmission.arguments.nsteps,
					"-m", jobsubmission.arguments.phase,
					// input files
					"-f", jobsubmission.jobdir+"/data.mzxml",
					"-s", jobsubmission.jobdir+"/smiles.sd",
					// output file
					"-o", jobsubmission.jobdir+"/results.db"
			};
	        sd.setArguments(args);
	        sd.setStderr(GAT.createFile("stderr.txt"));
            sd.setStdout(GAT.createFile("stdout.txt"));
            Preferences preferences = new Preferences();
            preferences.put("resourcebroker.adaptor.name", "local");
            wjd.add(new JobDescription(sd), new URI("any://localhost"),
                    preferences);

	        return wjd;
		} else {
			throw new Exception("Unknown job type: '" + jobsubmission.jobtype+"'");
		}
	}
}
