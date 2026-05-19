using PhysicsSim.Core;
using Xunit;

namespace PhysicsSim.Core.Tests;

public class AstronomyConstantsTests
{
    [Fact]
    public void AstronomicalUnit_IsAbout150MillionKm()
    {
        Assert.Equal(1.495978707e11, AstronomyConstants.AstronomicalUnit, 1);
    }

    [Fact]
    public void PlanetMasses_ArePositiveAndOrdered()
    {
        Assert.True(AstronomyConstants.MercuryMass > 0);
        Assert.True(AstronomyConstants.VenusMass > AstronomyConstants.MercuryMass);
        Assert.True(AstronomyConstants.EarthMass > AstronomyConstants.VenusMass);
        Assert.True(AstronomyConstants.JupiterMass > AstronomyConstants.EarthMass);
        Assert.True(AstronomyConstants.SolarMass > AstronomyConstants.JupiterMass);
    }

    [Fact]
    public void PlanetRadii_ArePositive()
    {
        Assert.True(AstronomyConstants.MercuryRadius > 0);
        Assert.True(AstronomyConstants.VenusRadius > 0);
        Assert.True(AstronomyConstants.EarthRadius > 0);
        Assert.True(AstronomyConstants.MarsRadius > 0);
        Assert.True(AstronomyConstants.JupiterMeanRadius > 0);
        Assert.True(AstronomyConstants.SaturnMeanRadius > 0);
    }

    [Fact]
    public void JupiterLowFlyby_IsTwiceRadius()
    {
        Assert.Equal(2.0 * AstronomyConstants.JupiterMeanRadius,
                     AstronomyConstants.JupiterLowFlybyDistance, 1);
    }

    [Fact]
    public void SemiMajorAxes_InCorrectRange()
    {
        Assert.True(AstronomyConstants.JupiterSemiMajorAxis > 4 * AstronomyConstants.AstronomicalUnit);
        Assert.True(AstronomyConstants.JupiterSemiMajorAxis < 6 * AstronomyConstants.AstronomicalUnit);
        Assert.True(AstronomyConstants.SaturnSemiMajorAxis > AstronomyConstants.JupiterSemiMajorAxis);
    }

    [Fact]
    public void GravitationalParameters_ArePositiveAndOrdered()
    {
        Assert.True(AstronomyConstants.MercuryGM > 0);
        Assert.True(AstronomyConstants.VenusGM > AstronomyConstants.MercuryGM);
        Assert.True(AstronomyConstants.EarthGM > AstronomyConstants.VenusGM);
        Assert.True(AstronomyConstants.JupiterGM > AstronomyConstants.EarthGM);
        Assert.True(AstronomyConstants.SolarGM > AstronomyConstants.JupiterGM);
    }

    [Fact]
    public void Mass_ConsistentWithGM()
    {
        double relTol = 1e-4;
        Assert.Equal(AstronomyConstants.JupiterGM,
            AstronomyConstants.G * AstronomyConstants.JupiterMass, AstronomyConstants.JupiterGM * relTol);
        Assert.Equal(AstronomyConstants.SolarGM,
            AstronomyConstants.G * AstronomyConstants.SolarMass, AstronomyConstants.SolarGM * relTol);
    }
}
