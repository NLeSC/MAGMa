package nl.esciencecenter.magma.gat;

import static com.yammer.dropwizard.testing.JsonHelpers.fromJson;
import static com.yammer.dropwizard.testing.JsonHelpers.jsonFixture;
import static org.fest.assertions.api.Assertions.assertThat;

import java.io.IOException;

import org.junit.Test;

import com.google.common.collect.ImmutableMap;

public class GATConfigurationTest {

    @Test
    public void testGATConfigurationBrokerPrefs() {
        ImmutableMap<String, Object> prefs = ImmutableMap.of(
                "localq.max.concurrent.jobs", (Object) 1);

        GATConfiguration g = new GATConfiguration("localq://localhost", prefs);

        assertThat(g.getBrokerURI()).isEqualTo("localq://localhost");
        assertThat(g.getPreferences()).isEqualTo(prefs);
    }

    @Test
    public void testSetGetBrokerURI() {
        GATConfiguration g = new GATConfiguration();
        g.setBrokerURI("localq://localhost");
        assertThat(g.getBrokerURI()).isEqualTo("localq://localhost");
    }

    @Test
    public void testGetPreferences() {
        GATConfiguration g = new GATConfiguration();
        ImmutableMap<String, Object> expected = ImmutableMap.of();
        assertThat(g.getPreferences()).isEqualTo(expected);
    }

    @Test
    public void testSetPreferences() {
        GATConfiguration g = new GATConfiguration();
        ImmutableMap<String, Object> prefs = ImmutableMap.of("mykey",
                (Object) "myval");
        g.setPreferences(prefs);
        assertThat(g.getPreferences()).isEqualTo(prefs);
    }

    @Test
    public void deserializesFromJSON() throws IOException {
        GATConfiguration actual = fromJson(jsonFixture("fixtures/gat.json"),
                GATConfiguration.class);

        assertThat(actual.getBrokerURI()).isEqualTo("localq://localhost");
        ImmutableMap<String, Object> prefs = ImmutableMap.of(
                "localq.max.concurrent.jobs", (Object) 1);
        assertThat(actual.getPreferences()).isEqualTo(prefs);
    }

    @Test
    public void hasAWorkingEqualsMethod() throws Exception {
        ImmutableMap<String, Object> prefs = ImmutableMap.of(
                "localq.max.concurrent.jobs", (Object) 1);

        GATConfiguration g = new GATConfiguration("localq://localhost", prefs);

        assertThat(g.equals(g)).isTrue();

        assertThat(g.equals(new GATConfiguration("localq://localhost", prefs)))
                .isTrue();

        assertThat(g.equals(null)).isFalse();

        assertThat(g.equals("string")).isFalse();

        assertThat(g.equals(new GATConfiguration("ssh://example.com", prefs)))
                .isFalse();

        ImmutableMap<String, Object> prefs2 = ImmutableMap.of(
                "localq.max.concurrent.jobs", (Object) 4);
        assertThat(g.equals(new GATConfiguration("localq://localhost", prefs2)))
                .isFalse();
    }
}
