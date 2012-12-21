package nl.esciencecenter.magma.api;

import com.google.common.base.Objects;

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

	@Override
	public boolean equals(Object obj) {
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		JobSubmitResponse other = (JobSubmitResponse) obj;
		return Objects.equal(this.jobid, other.jobid);
	}
}
