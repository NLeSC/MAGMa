package org.nlesc.magma.entities;

import javax.xml.bind.annotation.XmlRootElement;

@XmlRootElement
public class JobSubmitRequest {
    public String jobdir;
    public String jobtype;
    // wanted to use Map<String,String> as arguments type but jackson/jersey can't work with it
    // so switched to explicit class making it impossible to reuse for jobtypes with different arguments
    public MzxmlArguments arguments;

    public JobSubmitRequest(String jobdir, String jobtype) {
        this.jobdir = jobdir;
        this.jobtype = jobtype;
        this.arguments = new MzxmlArguments();
    }

    public JobSubmitRequest(String jobdir, String jobtype, MzxmlArguments arguments) {
        this.jobdir = jobdir;
        this.jobtype = jobtype;
        this.arguments = arguments;
    }

    public JobSubmitRequest() {}
}
