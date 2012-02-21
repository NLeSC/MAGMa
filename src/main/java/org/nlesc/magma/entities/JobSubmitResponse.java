package org.nlesc.magma.entities;

import javax.xml.bind.annotation.XmlRootElement;

/**
 * Returned when job is submitted successfully.
 *
 * @author Stefan Verhoeven <s.verhoeven@esciencecenter.nl>
 */
@XmlRootElement
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
