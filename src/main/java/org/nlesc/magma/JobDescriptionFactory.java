package org.nlesc.magma;

import org.gridlab.gat.GAT;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.SoftwareDescription;
import org.nlesc.magma.entities.JobSubmitRequest;

public class JobDescriptionFactory {

    /**
     * Based on jobsubmission.jobtype returns a job description.
     *
     * @param jobsubmission Job submission request
     * @return JobDescription
     * @throws Exception when jobtype is not implemented
     */
    public JobDescription getJobDescription(JobSubmitRequest jobsubmission) throws Exception {
        if (jobsubmission.jobtype.equals("sleep")) {
                SoftwareDescription sd = new SoftwareDescription();
                // simulate job by sleeping a while
                sd.setExecutable("/bin/sleep");
                String[] args = {"30"};
                sd.setArguments(args);

                JobDescription jd = new JobDescription(sd);
                return jd;
        } else if (jobsubmission.jobtype.equals("mzxmllocal")) {

            SoftwareDescription sd = new SoftwareDescription();
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
                    "-f", jobsubmission.jobdir + "/data.mzxml",
                    "-s", jobsubmission.jobdir + "/smiles.sd",
                    // output file
                    "-o", jobsubmission.jobdir + "/results.db"
            };
            sd.setArguments(args);

            sd.setStdout(GAT.createFile(jobsubmission.jobdir + "/stdout.txt"));
            sd.setStderr(GAT.createFile(jobsubmission.jobdir + "/stderr.txt"));

            JobDescription jd = new JobDescription(sd);
            return jd;
        } else if (jobsubmission.jobtype.equals("mzxmlremote")) {
            SoftwareDescription sd = new SoftwareDescription();

            sd.setExecutable("/bin/sh");
            String [] args = {
                    "mscore_mzxml.sh", // expects mscore_mzxml to be in path
                    // "echo",
                    "-p", jobsubmission.arguments.precision,
                    "-c", jobsubmission.arguments.mscutoff,
                    "-d", jobsubmission.arguments.msmscutoff,
                    "-i", jobsubmission.arguments.ionisation,
                    "-n", jobsubmission.arguments.nsteps,
                    "-m", jobsubmission.arguments.phase,
                    // input files
                    "-f", "data.mzxml",
                    "-s", "smiles.sd",
                    // output file
                    "-o", "results.db"
            };
            sd.setArguments(args);

            sd.setStdout(GAT.createFile(jobsubmission.jobdir + "/stdout.txt"));
            sd.setStderr(GAT.createFile(jobsubmission.jobdir + "/stderr.txt"));

            sd.addPreStagedFile(GAT.createFile("Magma-1.1.tar.gz"));
            sd.addPreStagedFile(GAT.createFile("mscore_mzxml.sh"));
            sd.addPreStagedFile(GAT.createFile(jobsubmission.jobdir + "/data.mzxml"));
            sd.addPreStagedFile(GAT.createFile(jobsubmission.jobdir + "/smiles.sd"));

            sd.addPostStagedFile(GAT.createFile("results.db"),
                    GAT.createFile(jobsubmission.jobdir + "/results.db"));

            JobDescription jd = new JobDescription(sd);
            return jd;
        } else {
            throw new Exception("Unknown job type: '" + jobsubmission.jobtype + "'");
        }
    }
}
