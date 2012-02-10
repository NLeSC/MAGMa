package org.nlesc.magma;


import java.util.UUID;

import static org.junit.Assert.*;

import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.WrapperJobDescription;
import org.junit.Test;
import org.nlesc.magma.entities.JobSubmitRequest;
import org.nlesc.magma.entities.MzxmlArguments;

import junit.framework.TestCase;

public class JobDescriptionFactoryTest extends TestCase {

	@Test
	public void testGetJobDescriptionSleep() {
		JobDescriptionFactory fact = new JobDescriptionFactory();

		JobSubmitRequest jobsubmission = new JobSubmitRequest("/somepath", "sleep");
		try {
			JobDescription jd = fact.getJobDescription(jobsubmission);
			assertEquals(jd.getSoftwareDescription().getExecutable(), "/bin/sleep");
			assertSame(jd.getSoftwareDescription().getArguments()[0], "30");
		} catch (Exception e) {
			fail("Sleep job type must exist");
		}
	}

	@Test
	public void testGetJobDescriptionUnknown() {
		JobDescriptionFactory fact = new JobDescriptionFactory();
		JobSubmitRequest jobsubmission = new JobSubmitRequest("/somepath", "holdbreath");
		try {
			fact.getJobDescription(jobsubmission);
			fail("Holdbreath job type must not exist");
		} catch (Exception e) {
			assertEquals(e.getMessage(), "Unknown job type: 'holdbreath'");
		}
	}

	@Test
	public void testGetJobDescriptionMzxmlLocal() {
		JobDescriptionFactory fact = new JobDescriptionFactory();

		// When user submits a magma job the website will put
		// the user supplied mzxml file, smiles file and config file into a jobdir.
		// then do a rest call to this service
		// which should if broker is local
		// 1. run mscore_mzxml with parameters in config file and files in jobdir
		// 2. run sygma with parameters in config file
		// 3. write completed flag to jobdir
		// if broker is glite then
		// 1. sync jobdir to srm
		// 2. submit wrapper job and inside wrapper job do
		// 2.1 sync jobdir from srm to localdir
		// 2.2 fetch program tarball from srm
		// 2.3 unpack program tarball
		// 2.3.1 run progs same as for local broker
		// 2.4 sync localdir to srm
		// 3. sync jobdir from srm
		// 4. clean srm
		// 5. write completed flag to jobdir

		String jobdir = System.getProperty("java.io.tmpdir")+"/"+UUID.randomUUID().toString();

		MzxmlArguments arguments = new MzxmlArguments();
		// arguments for mzxml
		arguments.precision= "0.01";
		arguments.mscutoff = "2e5";
		arguments.msmscutoff = "0.1";
		arguments.ionisation = "pos";
		// arguments for sygma
		arguments.nsteps = "2";
		arguments.phase = "12";

		JobSubmitRequest jobsubmission = new JobSubmitRequest(jobdir, "mzxmllocal", arguments);
		try {
			JobDescription jd = fact.getJobDescription(jobsubmission);

			assertEquals(jd.getSoftwareDescription().getExecutable(), "/usr/bin/env");
			String[] jobargs = {
					"mscore_mzxml",
					"-p", arguments.precision,
					"-c", arguments.mscutoff,
					"-d", arguments.msmscutoff,
					"-i", arguments.ionisation,
					"-n", arguments.nsteps,
					"-m", arguments.phase,
					"-f", jobdir+"/data.mzxml",
					"-s", jobdir+"/smiles.sd",
					"-o", jobdir+"/results.db"
			};
			assertArrayEquals(jd.getSoftwareDescription().getArguments(), jobargs);
			assertEquals(jd.getSoftwareDescription().getStdout().getAbsolutePath(), jobdir+"/stdout.txt");
			assertEquals(jd.getSoftwareDescription().getStderr().getAbsolutePath(), jobdir+"/stderr.txt");
		} catch (Exception e) {
			fail("Mzxml local job type must exist");
		}
	}

	@Test
	public void testGetJobDescriptionMzxmlRemote() {
		JobDescriptionFactory fact = new JobDescriptionFactory();

		// When user submits a magma job the website will put
		// the user supplied mzxml file, smiles file and config file into a jobdir.
		// then do a rest call to this service
		// which should if broker is local
		// 1. run mscore_mzxml with parameters in config file and files in jobdir
		// 2. run sygma with parameters in config file
		// 3. write completed flag to jobdir
		// if broker is remote then
		// 1. sync jobdir to srm
		// 2. submit wrapper job and inside wrapper job do
		// 2.1 sync jobdir from srm to localdir
		// 2.2 fetch program tarball from srm
		// 2.3 unpack program tarball
		// 2.3.1 run progs same as for local broker
		// 2.4 sync localdir to srm
		// 3. sync jobdir from srm
		// 4. clean srm
		// 5. write completed flag to jobdir

		String jobdir = System.getProperty("java.io.tmpdir")+"/"+UUID.randomUUID().toString();

		MzxmlArguments arguments = new MzxmlArguments();
		// arguments for mzxml
		arguments.precision= "0.01";
		arguments.mscutoff = "2e5";
		arguments.msmscutoff = "0.1";
		arguments.ionisation = "pos";
		// arguments for sygma
		arguments.nsteps = "2";
		arguments.phase = "12";

		JobSubmitRequest jobsubmission = new JobSubmitRequest(jobdir, "mzxmlremote", arguments);
		try {
			WrapperJobDescription wrappedjd = (WrapperJobDescription) fact.getJobDescription(jobsubmission);
			JobDescription jd = wrappedjd.getJobInfos().get(0).getJobDescription();

			assertEquals(jd.getSoftwareDescription().getExecutable(), "/usr/bin/env");
			String[] jobargs = {
					"mscore_mzxml",
					"-p", arguments.precision,
					"-c", arguments.mscutoff,
					"-d", arguments.msmscutoff,
					"-i", arguments.ionisation,
					"-n", arguments.nsteps,
					"-m", arguments.phase,
					"-f", jobdir+"/data.mzxml",
					"-s", jobdir+"/smiles.sd",
					"-o", jobdir+"/results.db"
			};
			assertArrayEquals(jd.getSoftwareDescription().getArguments(), jobargs);
		} catch (Exception e) {
			fail("Mzxml local job type must exist");
		}
	}
}
