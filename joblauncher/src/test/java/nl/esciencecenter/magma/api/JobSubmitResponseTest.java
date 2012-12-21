package nl.esciencecenter.magma.api;

import static com.yammer.dropwizard.testing.JsonHelpers.asJson;
import static com.yammer.dropwizard.testing.JsonHelpers.fromJson;
import static com.yammer.dropwizard.testing.JsonHelpers.jsonFixture;
import static org.hamcrest.Matchers.equalTo;
import static org.hamcrest.Matchers.is;
import static org.junit.Assert.*;
import static org.fest.assertions.api.Assertions.assertThat;

import java.io.IOException;

import org.junit.Test;

public class JobSubmitResponseTest {

	private JobSubmitResponse sampleResponse() {
		return new JobSubmitResponse("1234");
	}

	@Test
	public void testJobSubmitResponseString() {
		assertEquals("1234", sampleResponse().jobid);
	}

	@Test
	public void testJobSubmitResponse() {
		JobSubmitResponse r = new JobSubmitResponse();
		assertNull(r.jobid);
	}

	@Test
	public void serializesToJSON() throws IOException {
		assertThat("a JobSubmitResponse can be serialized to JSON",
				asJson(sampleResponse()),
				is(equalTo(jsonFixture("fixtures/response.json"))));
	}

	@Test
	public void deserializesFromJSON() throws IOException {
		assertThat(
				"a JobSubmitResponse can be deserialized from JSON",
				fromJson(jsonFixture("fixtures/response.json"),
						JobSubmitResponse.class), is(sampleResponse()));
	}

	@Test
	public void hasAWorkingEqualsMethod() throws Exception {
		JobSubmitResponse response = sampleResponse();
		assertThat(response.equals(response)).isTrue();

		assertThat(response.equals(new JobSubmitResponse("1234"))).isTrue();

		assertThat(response.equals(null)).isFalse();

		assertThat(response.equals("string")).isFalse();

		assertThat(response.equals(new JobSubmitResponse("u1"))).isFalse();
	}

}
