using PhysicsSim.Core;
using PhysicsSim.Core.Ode;
using Xunit;

namespace PhysicsSim.Core.Tests;

public class NBodySimulatorTests
{
    private static IntegrationSettings FastSettings => new()
    {
        AbsTol = 1e-8,
        RelTol = 1e-8,
        InitialStep = 0.01,
        MinStep = 1e-10,
        MaxStep = 0.5,
    };

    [Fact]
    public void CircularOrbit_EnergyConserved()
    {
        double G = 1.0;
        double M = 1000.0;
        double m = 1.0;
        double r = 10.0;
        double v = Math.Sqrt(G * M / r);

        var bodies = new List<BodyState>
        {
            new("Star", M, Vector3d.Zero, Vector3d.Zero),
            new("Planet", m, new Vector3d(r, 0, 0), new Vector3d(0, v, 0)),
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: true);
        double period = 2 * Math.PI * r / v;

        var result = NBodySimulator.Simulate(
            system, t0: 0, t1: period, outputDt: period / 20,
            settings: FastSettings);

        Assert.True(result.SampleCount > 10);
        Assert.Null(result.Collision);

        double e0 = result.TotalEnergies![0];
        double eFinal = result.TotalEnergies[^1];
        double relError = Math.Abs((eFinal - e0) / Math.Abs(e0));
        Assert.True(relError < 1e-6, $"Energy drift {relError:E2} too large");
    }

    [Fact]
    public void CircularOrbit_MomentumConserved()
    {
        double G = 1.0;
        double M = 1000.0;
        double m = 1.0;
        double r = 10.0;
        double v = Math.Sqrt(G * M / r);

        var bodies = new List<BodyState>
        {
            new("Star", M, Vector3d.Zero, Vector3d.Zero),
            new("Planet", m, new Vector3d(r, 0, 0), new Vector3d(0, v, 0)),
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: true);
        double period = 2 * Math.PI * r / v;

        var result = NBodySimulator.Simulate(
            system, t0: 0, t1: period, outputDt: period / 20,
            settings: FastSettings);

        var p0 = result.TotalMomenta![0];
        var pFinal = result.TotalMomenta[^1];
        double drift = (pFinal - p0).Length();
        Assert.True(drift < 1e-8, $"Momentum drift {drift:E2} too large");
    }

    [Fact]
    public void SimulationResult_ContainsExpectedData()
    {
        double G = 1.0;
        var bodies = new List<BodyState>
        {
            new("A", 100.0, Vector3d.Zero, Vector3d.Zero),
            new("B", 1.0, new Vector3d(5, 0, 0), new Vector3d(0, 4, 0)),
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: false);
        var result = NBodySimulator.Simulate(
            system, 0, 2, 0.5, FastSettings);

        Assert.Equal(2, result.BodyCount);
        Assert.Equal("A", result.BodyNames[0]);
        Assert.Equal("B", result.BodyNames[1]);
        Assert.True(result.SampleCount >= 4);
        Assert.Equal(result.SampleCount, result.Positions.Length);
        Assert.Equal(result.SampleCount, result.Velocities.Length);
        Assert.NotNull(result.TotalEnergies);
        Assert.NotNull(result.TotalMomenta);
    }

    [Fact]
    public void CollisionDetection_DetectsCollision()
    {
        double G = 1.0;
        var bodies = new List<BodyState>
        {
            new("A", 100.0, new Vector3d(0, 0, 0), new Vector3d(0, 0, 0)),
            new("B", 1.0, new Vector3d(2, 0, 0), new Vector3d(-5, 0, 0)),
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: false);
        var radii = new List<double> { 0.5, 0.5 };

        var result = NBodySimulator.Simulate(
            system, 0, 2, 0.05, FastSettings,
            collisionRadii: radii);

        Assert.NotNull(result.Collision);
        Assert.True(result.TerminatedEarly);
        Assert.Equal("A", result.Collision.BodyAName);
        Assert.Equal("B", result.Collision.BodyBName);
    }

    [Fact]
    public void ThreeBodySystem_Integrates()
    {
        double G = 1.0;
        var bodies = new List<BodyState>
        {
            new("A", 100.0, new Vector3d(0, 0, 0), new Vector3d(0, 0, 0)),
            new("B", 10.0, new Vector3d(5, 0, 0), new Vector3d(0, 4, 0)),
            new("C", 1.0, new Vector3d(0, 8, 0), new Vector3d(-3, 0, 0)),
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: true);
        var result = NBodySimulator.Simulate(
            system, 0, 2, 0.5, FastSettings);

        Assert.Equal(3, result.BodyCount);
        Assert.True(result.SampleCount >= 4);

        double e0 = result.TotalEnergies![0];
        double eFinal = result.TotalEnergies[^1];
        double relError = Math.Abs((eFinal - e0) / Math.Abs(e0));
        Assert.True(relError < 1e-4, $"Three-body energy drift {relError:E2}");
    }

    [Fact]
    public void StopCondition_TerminatesSimulation()
    {
        double G = 1.0;
        var bodies = new List<BodyState>
        {
            new("A", 100.0, Vector3d.Zero, Vector3d.Zero),
            new("B", 1.0, new Vector3d(5, 0, 0), new Vector3d(0, 20, 0)),
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: false);
        var result = NBodySimulator.Simulate(
            system, 0, 10, 0.1, FastSettings,
            stopCondition: (t, state) =>
            {
                double x = state[6];
                double y = state[7];
                double dist = Math.Sqrt(x * x + y * y);
                return dist > 15 ? "Escaped" : null;
            });

        Assert.NotNull(result.TerminationReason);
        Assert.True(result.TerminatedEarly);
    }
}
