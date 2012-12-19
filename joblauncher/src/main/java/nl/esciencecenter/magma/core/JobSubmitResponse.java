package nl.esciencecenter.magma.core;

/**
 * Returned when job is submitted successfully.
 *
 * @author Stefan Verhoeven <s.verhoeven@esciencecenter.nl>
 */
public class JobSubmitResponse {
    public String jobid;

    public JobSubmitResponse(String jobid) {
        this.jobid = jobid;
    }

    public JobSubmitResponse() {
        super();
        jobid = null;
    }
}
