package org.nlesc.magma;


import static org.junit.Assert.assertArrayEquals;

import java.util.HashMap;
import java.util.UUID;

import junit.framework.TestCase;

import org.gridlab.gat.GAT;
import org.gridlab.gat.io.File;
import org.gridlab.gat.resources.JobDescription;
import org.gridlab.gat.resources.SoftwareDescription;
import org.junit.Test;
import org.nlesc.magma.entities.JobSubmitRequest;
import org.nlesc.magma.entities.MzxmlArguments;

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
        // the user supplied mzxml file and smiles file into a jobdir.
        // then do a rest call to this service
        // which should if broker is local
        // 1. run mscore_mzxml with parameters in config file and files in jobdir
        // 2. run sygma with parameters from rest arguments
        // 3. write completed flag to jobdir

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
        // which if broker is remote then
        // 1. Pre stage input files
        // 2. Pre stage program tarball
        // 3. Pre stage shell script
        // 4. Run shell script
        // 4.1 unpack program tarball
        // 4.2 run progs same as for local broker
        // 5. Post stage result file
        // 6. write completed flag to jobdir
        //
        // shell script contains:
        // tar -zxf Magma-1.1.tar.gz
        // Magma-1.1/mscore_mzxml.sh $@

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
            JobDescription jd = fact.getJobDescription(jobsubmission);

            SoftwareDescription sd = jd.getSoftwareDescription();
            assertEquals(sd.getExecutable(), "/bin/sh");
            String[] jobargs = {
                    "mscore_mzxml.sh",
                    "-p", arguments.precision,
                    "-c", arguments.mscutoff,
                    "-d", arguments.msmscutoff,
                    "-i", arguments.ionisation,
                    "-n", arguments.nsteps,
                    "-m", arguments.phase,
                    "-f", "data.mzxml",
                    "-s", "smiles.sd",
                    "-o", "results.db"
            };
            assertArrayEquals(sd.getArguments(), jobargs);
            assertEquals(sd.getStdout().getAbsolutePath(), jobdir+"/stdout.txt");
            assertEquals(sd.getStderr().getAbsolutePath(), jobdir+"/stderr.txt");
            /**
             * assert pre and post staged files
             * Pre:
             * {
             * "Magma-1.1.tar.gz": null,
             * "mscore_mzxml.sh": null,
             * jobdir+"data.mzxml": null,
             * jobdir+"smiles.sd": null
             * }
             *
             * Post:
             * { "results.db": jobdir+"/results.db" }
             */
            HashMap<File, File> expectedPreStaged = new HashMap<File, File>();
            expectedPreStaged.put(GAT.createFile("Magma-1.1.tar.gz"), null);
            expectedPreStaged.put(GAT.createFile("mscore_mzxml.sh"), null);
            expectedPreStaged.put(GAT.createFile(jobdir+"/data.mzxml"), null);
            expectedPreStaged.put(GAT.createFile(jobdir+"/smiles.sd"), null);
            assertEquals(expectedPreStaged, sd.getPreStaged());

            HashMap<File, File> expectedPostStaged = new HashMap<File, File>();
            expectedPostStaged.put(GAT.createFile("results.db"), GAT.createFile(jobdir+"/results.db"));
            assertEquals(expectedPostStaged, sd.getPostStaged());
        } catch (Exception e) {
            fail("Mzxml remote job type must exist");
        }
    }
}
