namespace PhysicsSim.Core;

public static class AstronomyConstants
{
    public const double G = 6.67430e-11;
    public const double AstronomicalUnit = 1.495978707e11;

    // GM (гравитационные параметры, м³/с²) — точнее, чем G*M по отдельности.
    // Источник: IAU 2015 / JPL DE440.
    public const double SolarGM        = 1.32712440041279419e20;
    public const double MercuryGM      = 2.2031868551e13;
    public const double VenusGM        = 3.24858592000e14;
    public const double EarthGM        = 3.98600435507e14;
    public const double MarsGM         = 4.2828375816e13;
    public const double JupiterGM      = 1.26712764100000e17;
    public const double SaturnGM       = 3.79395259000000e16;
    public const double UranusGM       = 5.79393658000000e15;
    public const double NeptuneGM      = 6.83509920000000e15;

    // Массы (кг) — для справки и совместимости (вычислены как GM/G).
    public const double SolarMass = SolarGM / G;
    public const double MercuryMass = MercuryGM / G;
    public const double VenusMass = VenusGM / G;
    public const double EarthMass = EarthGM / G;
    public const double MarsMass = MarsGM / G;
    public const double JupiterMass = JupiterGM / G;
    public const double SaturnMass = SaturnGM / G;
    public const double UranusMass = UranusGM / G;
    public const double NeptuneMass = NeptuneGM / G;

    public const double SolarRadius = 6.9634e8;
    public const double MercuryRadius = 2.4397e6;
    public const double VenusRadius = 6.0518e6;
    public const double EarthRadius = 6.3710e6;
    public const double MarsRadius = 3.3895e6;
    public const double JupiterMeanRadius = 7.1492e7;
    public const double SaturnMeanRadius = 5.8232e7;
    public const double UranusRadius = 2.5362e7;
    public const double NeptuneRadius = 2.4622e7;
    public const double JupiterLowFlybyDistance = 2.0 * JupiterMeanRadius;

    public const double JupiterSemiMajorAxis = 5.2044 * AstronomicalUnit;
    public const double SaturnSemiMajorAxis = 9.5826 * AstronomicalUnit;
}
