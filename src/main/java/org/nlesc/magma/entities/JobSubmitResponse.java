package org.nlesc.magma.entities;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement
public class JobSubmitResponse {
	public String jobid;

	public JobSubmitResponse(String jobid) {
		this.jobid = jobid;
	}

	public JobSubmitResponse() {
		super();
		// TODO Auto-generated constructor stub
	}
}
