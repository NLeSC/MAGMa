package nl.nlesc.magma;

import nl.nlesc.magma.entities.JobSubmitRequest;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.io.File;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.SoftwareDescription;

public class JobDescriptionFactory {

    /**
     * Convert requested jobsubmission to JobDescription which can be submitted
     *
     * @param jobsubmission
     * @return JobDescription
     * @throws GATObjectCreationException
     */
    public JobDescription getJobDescription(
            JobSubmitRequest jobsubmission) throws GATObjectCreationException {

        SoftwareDescription sd = new SoftwareDescription();
        sd.setExecutable(jobsubmission.executable);
        sd.setArguments(jobsubmission.arguments);
        sd.setStderr(GAT.createFile(jobsubmission.jobdir + jobsubmission.stderr));
        sd.setStdout(GAT.createFile(jobsubmission.jobdir + jobsubmission.stdout));
        if (jobsubmission.time_max > 0) {
            sd.addAttribute("time.max", jobsubmission.time_max);
        }
        if (jobsubmission.memory_min > 0) {
            sd.addAttribute("memory.min", jobsubmission.memory_min);
        }
        if (jobsubmission.memory_max > 0) {
            sd.addAttribute("memory.max", jobsubmission.memory_max);
        }
        for (String prestage: jobsubmission.prestaged) {
            File prestagefile = GAT.createFile(prestage);
            if (!prestagefile.isAbsolute()) {
                prestagefile = GAT.createFile(jobsubmission.jobdir + prestage);
            }
            sd.addPreStagedFile(prestagefile);
        }
        for (String poststage: jobsubmission.poststaged) {
            File poststagefile = GAT.createFile(poststage);
            sd.addPostStagedFile(poststagefile, GAT.createFile(jobsubmission.jobdir + poststage));
        }

        return new JobDescription(sd);
    }
}
