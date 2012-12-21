package nl.esciencecenter.magma.api;

import java.net.URI;
import java.util.Arrays;

import javax.validation.constraints.NotNull;

import org.gridlab.gat.GAT;
import org.gridlab.gat.GATObjectCreationException;
import org.gridlab.gat.io.File;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.SoftwareDescription;

/**
 * Request which can be converted to JobDescription which can be submitted using
 * JavaGAT.
 *
 * @author Stefan Verhoeven <s.verhoeven@esciencecenter.nl>
 *
 */
public class JobSubmitRequest {
	/**
	 * Job directory where stderr/stdout/prestaged/poststaged file are relative
	 * to and where job.state file is written. Must end with '/'
	 */
	@NotNull
	public String jobdir;
	/**
	 * Path to executable on execution host
	 */
	@NotNull
	public String executable;
	/**
	 * File name to write standard error to
	 */
	public String stderr = "stderr.txt";
	/**
	 * File name to write standard out to
	 */
	public String stdout = "stdout.txt";
	/**
	 * Arguments passed to executable
	 */
	public String[] arguments = {};
	/**
	 * List of filenames to copy from job directory to work directory before
	 * executable is called. Work directory is created on the execution host.
	 * Can be relative to job directory or absolute paths.
	 */
	public String[] prestaged = {};
	/**
	 * List of filenames to copy from work directory to job directory after
	 * executable is called. Must be relative to job directory
	 */
	public String[] poststaged = {};
	/**
	 * The maximum walltime or cputime for a single execution of the executable.
	 * The units is in minutes.
	 */
	public long time_max = 0;
	/**
	 * minimal required memory in MB
	 */
	public int memory_min = 0;
	/**
	 * maximum required memory in MB
	 */
	public int memory_max = 0;
	/**
	 * Url where changes of state are PUT to.
	 */
	public URI status_callback_url;

	/**
	 * Constructor
	 *
	 * @param jobdir
	 * @param executable
	 * @param arguments
	 * @param prestaged
	 * @param poststaged
	 * @param stderr
	 * @param stdout
	 */
	public JobSubmitRequest(String jobdir, String executable,
			String[] arguments, String[] prestaged, String[] poststaged,
			String stderr, String stdout, URI status_callback_url) {
		super();
		this.jobdir = jobdir;
		this.executable = executable;
		this.arguments = arguments;
		this.prestaged = prestaged;
		this.poststaged = poststaged;
		this.stderr = stderr;
		this.stdout = stdout;
		this.status_callback_url = status_callback_url;
	}

	public URI getStatus_callback_url() {
		return status_callback_url;
	}

	/**
	 * JAXB needs this
	 */
	public JobSubmitRequest() {
		super();
	}

	/**
	 * Convert requested jobsubmission to JobDescription which can be submitted
	 *
	 * @return JobDescription
	 * @throws GATObjectCreationException
	 */
	public JobDescription toJobDescription() throws GATObjectCreationException {
		SoftwareDescription sd = new SoftwareDescription();
		sd.setExecutable(executable);
		sd.setArguments(arguments);
		sd.setStderr(GAT.createFile(jobdir + stderr));
		sd.setStdout(GAT.createFile(jobdir + stdout));
		if (time_max > 0) {
			sd.addAttribute("time.max", time_max);
		}
		if (memory_min > 0) {
			sd.addAttribute("memory.min", memory_min);
		}
		if (memory_max > 0) {
			sd.addAttribute("memory.max", memory_max);
		}
		for (String prestage : prestaged) {
			File prestagefile = GAT.createFile(prestage);
			if (!prestagefile.isAbsolute()) {
				prestagefile = GAT.createFile(jobdir + prestage);
			}
			sd.addPreStagedFile(prestagefile);
		}
		for (String poststage : poststaged) {
			File poststagefile = GAT.createFile(poststage);
			sd.addPostStagedFile(poststagefile,
					GAT.createFile(jobdir + poststage));
		}

		return new JobDescription(sd);
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		JobSubmitRequest other = (JobSubmitRequest) obj;
		if (!Arrays.equals(arguments, other.arguments))
			return false;
		if (executable == null) {
			if (other.executable != null)
				return false;
		} else if (!executable.equals(other.executable))
			return false;
		if (jobdir == null) {
			if (other.jobdir != null)
				return false;
		} else if (!jobdir.equals(other.jobdir))
			return false;
		if (memory_max != other.memory_max)
			return false;
		if (memory_min != other.memory_min)
			return false;
		if (!Arrays.equals(poststaged, other.poststaged))
			return false;
		if (!Arrays.equals(prestaged, other.prestaged))
			return false;
		if (status_callback_url == null) {
			if (other.status_callback_url != null)
				return false;
		} else if (!status_callback_url.equals(other.status_callback_url))
			return false;
		if (stderr == null) {
			if (other.stderr != null)
				return false;
		} else if (!stderr.equals(other.stderr))
			return false;
		if (stdout == null) {
			if (other.stdout != null)
				return false;
		} else if (!stdout.equals(other.stdout))
			return false;
		if (time_max != other.time_max)
			return false;
		return true;
	}

	@Override
	public String toString() {
		return "JobSubmitRequest [jobdir=" + jobdir + ", executable="
				+ executable + ", stderr=" + stderr + ", stdout=" + stdout
				+ ", arguments=" + Arrays.toString(arguments) + ", prestaged="
				+ Arrays.toString(prestaged) + ", poststaged="
				+ Arrays.toString(poststaged) + ", time_max=" + time_max
				+ ", memory_min=" + memory_min + ", memory_max=" + memory_max
				+ ", status_callback_url=" + status_callback_url + "]";
	}
}
