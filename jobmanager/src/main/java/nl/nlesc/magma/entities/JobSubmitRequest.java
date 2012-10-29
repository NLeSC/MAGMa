package nl.nlesc.magma.entities;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * Request which can be converted to JobDescription which can be submitted using
 * JavaGAT.
 *
 * @author Stefan Verhoeven <s.verhoeven@esciencecenter.nl>
 *
 */
@XmlRootElement
public class JobSubmitRequest {
    /**
     * Job directory where stderr/stdout/prestaged/poststaged file are relative
     * to and where job.state file is written.
     * Must end with '/'
     */
    public String jobdir;
    /**
     * Path to executable on execution host
     */
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
     * The maximum walltime or cputime for a single execution of the
     * executable. The units is in minutes.
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
     * Constructor
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
            String stderr, String stdout) {
        super();
        this.jobdir = jobdir;
        this.executable = executable;
        this.arguments = arguments;
        this.prestaged = prestaged;
        this.poststaged = poststaged;
        this.stderr = stderr;
        this.stdout = stdout;
    }

    /**
     * JAXB needs this
     */
    public JobSubmitRequest() {
        super();
    }
}
