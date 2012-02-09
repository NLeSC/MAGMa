package org.nlesc.magma.entities;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement
public class JobSubmitRequest {
	public String jobdir;
	public String jobtype; // mzxml or metabolize

	public JobSubmitRequest(String jobdir, String jobtype) {
		this.jobdir = jobdir;
		this.jobtype = jobtype;
	}

	public JobSubmitRequest() {
		super();
		// TODO Auto-generated constructor stub
	}
}
