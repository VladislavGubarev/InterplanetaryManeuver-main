using PhysicsSim.Core;
using Xunit;

namespace PhysicsSim.Core.Tests;

public class FlybyAnalysisTests
{
    [Fact]
    public void SphereOfInfluence_KnownValues()
    {
        double planetMass = AstronomyConstants.JupiterMass;
        double centralMass = AstronomyConstants.SolarMass;
        double semiMajorAxis = AstronomyConstants.JupiterSemiMajorAxis;

        double soi = FlybyAnalysis.ComputeSphereOfInfluenceRadius(planetMass, centralMass, semiMajorAxis);

        Assert.True(soi > 4e10, $"Jupiter SOI should be > 4e10 m, got {soi:E2}");
        Assert.True(soi < 6e10, $"Jupiter SOI should be < 6e10 m, got {soi:E2}");
    }

    [Fact]
    public void SphereOfInfluence_ScalesWithSemiMajorAxis()
    {
        double r1 = FlybyAnalysis.ComputeSphereOfInfluenceRadius(1e27, 2e30, 5e11);
        double r2 = FlybyAnalysis.ComputeSphereOfInfluenceRadius(1e27, 2e30, 10e11);

        Assert.True(r2 > r1, "SOI should increase with semi-major axis");
    }

    [Fact]
    public void SphereOfInfluence_ScalesWithPlanetMass()
    {
        double r1 = FlybyAnalysis.ComputeSphereOfInfluenceRadius(1e26, 2e30, 5e11);
        double r2 = FlybyAnalysis.ComputeSphereOfInfluenceRadius(1e27, 2e30, 5e11);

        Assert.True(r2 > r1, "SOI should increase with planet mass");
    }

    [Fact]
    public void Compute_EmptyResult_Throws()
    {
        var emptyResult = new SimulationResult
        {
            Times = Array.Empty<double>(),
            BodyNames = Array.Empty<string>(),
            Positions = Array.Empty<Vector3d[]>(),
            Velocities = Array.Empty<Vector3d[]>(),
        };

        Assert.Throws<ArgumentException>(() =>
            FlybyAnalysis.Compute(emptyResult, 0, 1, 2, 1e10));
    }

    [Fact]
    public void Compute_InvalidIndices_Throw()
    {
        var result = CreateMinimalResult(3, 10);

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            FlybyAnalysis.Compute(result, 0, 1, -1, 1e10));

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            FlybyAnalysis.Compute(result, 0, -1, 2, 1e10));

        Assert.Throws<ArgumentOutOfRangeException>(() =>
            FlybyAnalysis.Compute(result, -1, 1, 2, 1e10));
    }

    [Fact]
    public void Compute_BasicMetrics_AreReasonable()
    {
        int nBodies = 3;
        int nSamples = 100;

        double soiRadius = 50.0;
        var times = new double[nSamples];
        var positions = new Vector3d[nSamples][];
        var velocities = new Vector3d[nSamples][];

        for (int i = 0; i < nSamples; i++)
        {
            double t = i * 0.1;
            times[i] = t;

            double scX = 80.0 - 1.6 * i;
            double scY = 10.0 * Math.Sin(i * 0.1);

            positions[i] = new Vector3d[nBodies];
            velocities[i] = new Vector3d[nBodies];

            positions[i][0] = new Vector3d(500, 0, 0); // Sun
            positions[i][1] = new Vector3d(0, 0, 0);   // Jupiter
            positions[i][2] = new Vector3d(scX, scY, 0); // Spacecraft

            velocities[i][0] = new Vector3d(0, 0, 0);
            velocities[i][1] = new Vector3d(0, 13, 0);
            velocities[i][2] = new Vector3d(-16, 2 * Math.Cos(i * 0.1), 0);
        }

        var result = new SimulationResult
        {
            Times = times,
            BodyNames = new[] { "Sun", "Jupiter", "Spacecraft" },
            Positions = positions,
            Velocities = velocities,
        };

        var metrics = FlybyAnalysis.Compute(result, sunIndex: 0, planetIndex: 1, spacecraftIndex: 2, planetSoiRadius: soiRadius);

        Assert.True(metrics.MinDistanceToJupiter >= 0);
        Assert.True(metrics.InitialHeliocentricSpeed > 0);
        Assert.True(metrics.FinalHeliocentricSpeed > 0);
        Assert.True(metrics.InitialJupiterRelativeSpeed > 0);
    }

    [Fact]
    public void FlybyMetrics_IsFeasibleFlyby_RequiresAllConditions()
    {
        var feasible = new FlybyMetrics
        {
            InitialDistanceToJupiter = 1e9,
            FinalDistanceToJupiter = 1e9,
            JupiterSoiRadius = 5e10,
            MinDistanceToJupiter = 1e8,
            MinDistanceToSaturn = double.PositiveInfinity,
            ClosestApproachAltitudeToJupiter = 1e8 - AstronomyConstants.JupiterMeanRadius,
            InitialHeliocentricSpeed = 10000,
            FinalHeliocentricSpeed = 15000,
            InitialJupiterRelativeSpeed = 8000,
            FinalJupiterRelativeSpeed = 8000,
            DeltaVGainHeliocentric = 5000,
            EqualDistanceIndex = 50,
            EntryIndex = 10,
            ExitIndex = 90,
            ClosestApproachIndex = 40,
            JupiterCollisionIndex = -1,
        };

        Assert.True(feasible.IsFeasibleFlyby);
        Assert.True(feasible.HasSphereOfInfluenceCrossing);
        Assert.True(feasible.HasReturnToInitialDistance);
        Assert.False(feasible.HasJupiterCollision);
    }

    [Fact]
    public void FlybyMetrics_Collision_NotFeasible()
    {
        var collision = new FlybyMetrics
        {
            InitialDistanceToJupiter = 1e9,
            FinalDistanceToJupiter = 1e9,
            JupiterSoiRadius = 5e10,
            MinDistanceToJupiter = AstronomyConstants.JupiterMeanRadius * 0.5,
            MinDistanceToSaturn = double.PositiveInfinity,
            ClosestApproachAltitudeToJupiter = AstronomyConstants.JupiterMeanRadius * 0.5 - AstronomyConstants.JupiterMeanRadius,
            InitialHeliocentricSpeed = 10000,
            FinalHeliocentricSpeed = 15000,
            InitialJupiterRelativeSpeed = 8000,
            FinalJupiterRelativeSpeed = 8000,
            DeltaVGainHeliocentric = 5000,
            EqualDistanceIndex = 50,
            EntryIndex = 10,
            ExitIndex = 90,
            ClosestApproachIndex = 40,
            JupiterCollisionIndex = 35,
        };

        Assert.True(collision.HasJupiterCollision);
        Assert.False(collision.IsFeasibleFlyby);
    }

    [Fact]
    public void FlybyMetrics_DangerousLowFlyby()
    {
        var low = new FlybyMetrics
        {
            InitialDistanceToJupiter = 1e9,
            FinalDistanceToJupiter = 1e9,
            JupiterSoiRadius = 5e10,
            MinDistanceToJupiter = AstronomyConstants.JupiterMeanRadius * 1.5,
            MinDistanceToSaturn = double.PositiveInfinity,
            ClosestApproachAltitudeToJupiter = AstronomyConstants.JupiterMeanRadius * 0.5,
            InitialHeliocentricSpeed = 10000,
            FinalHeliocentricSpeed = 15000,
            InitialJupiterRelativeSpeed = 8000,
            FinalJupiterRelativeSpeed = 8000,
            DeltaVGainHeliocentric = 5000,
            EqualDistanceIndex = 50,
            EntryIndex = 10,
            ExitIndex = 90,
            ClosestApproachIndex = 40,
            JupiterCollisionIndex = -1,
        };

        Assert.True(low.HasDangerouslyLowJupiterFlyby);
        Assert.False(low.HasJupiterCollision);
    }

    private static SimulationResult CreateMinimalResult(int bodyCount, int sampleCount)
    {
        var times = new double[sampleCount];
        var positions = new Vector3d[sampleCount][];
        var velocities = new Vector3d[sampleCount][];
        var names = new string[bodyCount];

        for (int i = 0; i < bodyCount; i++)
            names[i] = $"Body{i}";

        for (int i = 0; i < sampleCount; i++)
        {
            times[i] = i;
            positions[i] = new Vector3d[bodyCount];
            velocities[i] = new Vector3d[bodyCount];
            for (int j = 0; j < bodyCount; j++)
            {
                positions[i][j] = new Vector3d(j * 10 + i, 0, 0);
                velocities[i][j] = new Vector3d(1, 0, 0);
            }
        }

        return new SimulationResult
        {
            Times = times,
            BodyNames = names,
            Positions = positions,
            Velocities = velocities,
        };
    }
}
