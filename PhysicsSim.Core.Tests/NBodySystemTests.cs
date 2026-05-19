using PhysicsSim.Core;
using Xunit;

namespace PhysicsSim.Core.Tests;

public class NBodySystemTests
{
    [Fact]
    public void Constructor_ThrowsForLessThanTwoBodies()
    {
        var single = new List<BodyState>
        {
            new("A", 1.0, Vector3d.Zero, Vector3d.Zero)
        };

        Assert.Throws<ArgumentException>(() => new NBodySystem(1.0, single));
    }

    [Fact]
    public void Dimension_IsSixTimesBodyCount()
    {
        var bodies = new List<BodyState>
        {
            new("A", 1.0, new Vector3d(0, 0, 0), new Vector3d(0, 0, 0)),
            new("B", 1.0, new Vector3d(1, 0, 0), new Vector3d(0, 1, 0)),
            new("C", 1.0, new Vector3d(0, 1, 0), new Vector3d(1, 0, 0)),
        };

        var system = new NBodySystem(1.0, bodies, toBarycentricFrame: false);
        Assert.Equal(18, system.Dimension);
        Assert.Equal(3, system.BodyCount);
    }

    [Fact]
    public void InitialState_StoresPositionsAndVelocities()
    {
        var bodies = new List<BodyState>
        {
            new("A", 10.0, new Vector3d(1, 2, 3), new Vector3d(4, 5, 6)),
            new("B", 20.0, new Vector3d(7, 8, 9), new Vector3d(10, 11, 12)),
        };

        var system = new NBodySystem(1.0, bodies, toBarycentricFrame: false);

        var a = system.ReadBodyFromState(system.InitialState, 0);
        Assert.Equal("A", a.Name);
        Assert.Equal(10.0, a.Mass);
        Assert.Equal(1, a.Position.X);
        Assert.Equal(2, a.Position.Y);
        Assert.Equal(3, a.Position.Z);

        var b = system.ReadBodyFromState(system.InitialState, 1);
        Assert.Equal("B", b.Name);
        Assert.Equal(20.0, b.Mass);
    }

    [Fact]
    public void BarycentricFrame_ShiftsMomentumToZero()
    {
        var bodies = new List<BodyState>
        {
            new("Sun", 1000.0, new Vector3d(0, 0, 0), new Vector3d(0, 0, 0)),
            new("Planet", 1.0, new Vector3d(10, 0, 0), new Vector3d(0, 10, 0)),
        };

        var system = new NBodySystem(1.0, bodies, toBarycentricFrame: true);
        var momentum = system.CalculateTotalMomentum(system.InitialState);

        Assert.True(momentum.Length() < 1e-10,
            $"Expected near-zero momentum, got {momentum.Length()}");
    }

    [Fact]
    public void ComputeDerivatives_PositionDerivativeEqualsVelocity()
    {
        var bodies = new List<BodyState>
        {
            new("A", 1.0, new Vector3d(0, 0, 0), new Vector3d(3, 4, 5)),
            new("B", 1.0, new Vector3d(10, 0, 0), new Vector3d(-1, 2, 0)),
        };

        var system = new NBodySystem(1.0, bodies, toBarycentricFrame: false, softeningSquared: 0);
        var dy = new double[system.Dimension];
        system.ComputeDerivatives(0, system.InitialState, dy);

        Assert.Equal(3, dy[0], 14);
        Assert.Equal(4, dy[1], 14);
        Assert.Equal(5, dy[2], 14);
    }

    [Fact]
    public void ComputeDerivatives_GravitationalAcceleration_TwoBodies()
    {
        double G = 1.0;
        double m = 1.0;
        double r = 1e6;

        var bodies = new List<BodyState>
        {
            new("A", m, new Vector3d(-r / 2, 0, 0), Vector3d.Zero),
            new("B", m, new Vector3d(r / 2, 0, 0), Vector3d.Zero),
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: false, softeningSquared: 0);
        var dy = new double[system.Dimension];
        system.ComputeDerivatives(0, system.InitialState, dy);

        double expectedAx = G * m / (r * r);
        Assert.Equal(expectedAx, dy[3], 12);
        Assert.Equal(0, dy[4], 14);
        Assert.Equal(0, dy[5], 14);

        Assert.Equal(-expectedAx, dy[9], 12);
    }

    [Fact]
    public void CalculateTotalEnergy_TwoBodySystem()
    {
        double G = 1.0;
        double m1 = 100.0, m2 = 1.0;
        var bodies = new List<BodyState>
        {
            new("A", m1, new Vector3d(0, 0, 0), new Vector3d(0, 0, 0)),
            new("B", m2, new Vector3d(5, 0, 0), new Vector3d(0, 4, 0)),
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: false);
        double energy = system.CalculateTotalEnergy(system.InitialState);

        double expectedKinetic = 0.5 * m2 * 16.0;
        double expectedPotential = -G * m1 * m2 / 5.0;
        Assert.Equal(expectedKinetic + expectedPotential, energy, 10);
    }

    [Fact]
    public void CalculateTotalMomentum_NoBarycentric()
    {
        var bodies = new List<BodyState>
        {
            new("A", 2.0, Vector3d.Zero, new Vector3d(3, 0, 0)),
            new("B", 5.0, new Vector3d(1, 0, 0), new Vector3d(0, 4, 0)),
        };

        var system = new NBodySystem(1.0, bodies, toBarycentricFrame: false);
        var p = system.CalculateTotalMomentum(system.InitialState);

        Assert.Equal(6, p.X, 10);
        Assert.Equal(20, p.Y, 10);
        Assert.Equal(0, p.Z, 10);
    }

    [Fact]
    public void GmConstructor_ProducesCorrectGMs()
    {
        var bodies = new List<BodyState>
        {
            new("Sun", 1e30, Vector3d.Zero, Vector3d.Zero),
            new("Planet", 1e24, new Vector3d(1e11, 0, 0), new Vector3d(0, 3e4, 0)),
        };
        var gms = new List<double> { 6.674e19, 6.674e13 };

        var system = new NBodySystem(bodies, gms, toBarycentricFrame: false);
        Assert.Equal(6.674e19, system.GMs[0], 1e10);
        Assert.Equal(6.674e13, system.GMs[1], 1e4);
        Assert.Equal(2, system.BodyCount);
    }

    [Fact]
    public void WriteBodyState_RoundTrip()
    {
        var bodies = new List<BodyState>
        {
            new("A", 1.0, Vector3d.Zero, Vector3d.Zero),
            new("B", 1.0, Vector3d.Zero, Vector3d.Zero),
        };

        var system = new NBodySystem(1.0, bodies, toBarycentricFrame: false);
        var pos = new Vector3d(7, 8, 9);
        var vel = new Vector3d(10, 11, 12);
        NBodySystem.WriteBodyState(system.InitialState, 1, pos, vel);

        var read = system.ReadBodyFromState(system.InitialState, 1);
        Assert.Equal(pos, read.Position);
        Assert.Equal(vel, read.Velocity);
    }

    [Fact]
    public void TensorPath_MatchesScalarPath()
    {
        double G = 6.674e-11;
        var bodies = new List<BodyState>
        {
            new("Sun", 1.989e30, new Vector3d(0, 0, 0), new Vector3d(0, 0, 0)),
            new("Jupiter", 1.898e27, new Vector3d(7.785e11, 0, 0), new Vector3d(0, 13070, 0)),
            new("Saturn", 5.683e26, new Vector3d(1.434e12, 0, 0), new Vector3d(0, 9690, 0)),
        };

        var scalar = new NBodySystem(G, bodies, toBarycentricFrame: false, useTensors: false);
        var tensor = new NBodySystem(G, bodies, toBarycentricFrame: false, useTensors: true);

        var dyScalar = new double[scalar.Dimension];
        var dyTensor = new double[tensor.Dimension];

        scalar.ComputeDerivatives(0, scalar.InitialState, dyScalar);
        tensor.ComputeDerivatives(0, tensor.InitialState, dyTensor);

        for (int k = 0; k < dyScalar.Length; k++)
        {
            Assert.Equal(dyScalar[k], dyTensor[k], 6);
        }
    }

    [Fact]
    public void TensorPath_EnergyConserved()
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

        var system = new NBodySystem(G, bodies, toBarycentricFrame: true, softeningSquared: 0, useTensors: true);
        double period = 2 * Math.PI * r / v;

        var result = NBodySimulator.Simulate(
            system, t0: 0, t1: period, outputDt: period / 20,
            settings: new PhysicsSim.Core.Ode.IntegrationSettings
            {
                AbsTol = 1e-8, RelTol = 1e-8, InitialStep = 0.01,
                MinStep = 1e-10, MaxStep = 0.5,
            });

        double e0 = result.TotalEnergies![0];
        double eFinal = result.TotalEnergies[^1];
        double relError = Math.Abs((eFinal - e0) / Math.Abs(e0));
        Assert.True(relError < 1e-6, $"Tensor energy drift {relError:E2}");
    }
}
