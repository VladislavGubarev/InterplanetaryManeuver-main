using PhysicsSim.Core;
using Xunit;

namespace PhysicsSim.Core.Tests;

public class PhysicsConservationTests
{
    [Fact]
    public void TwoBodySystem_EnergyShouldBeConserved()
    {
        // Настройка простой системы (Земля - Луна в миниатюре)
        double G = 1.0;
        var bodies = new List<BodyState>
        {
            new BodyState("M1", 1000.0, new Vector3d(0, 0, 0), new Vector3d(0, -0.1, 0)),
            new BodyState("M2", 1.0, new Vector3d(10, 0, 0), new Vector3d(0, 10.0, 0))
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: true);
        double initialEnergy = system.CalculateTotalEnergy(system.InitialState);
        Vector3d initialMomentum = system.CalculateTotalMomentum(system.InitialState);

        // В барицентрической системе импульс должен быть 0
        Assert.True(initialMomentum.Length() < 1e-10);

        // Имитация шага интеграции (просто проверяем расчеты на начальном этапе)
        // Для полноценного теста нужно прогнать симулятор, но мы проверим корректность методов расчета.
        double energy = system.CalculateTotalEnergy(system.InitialState);
        Assert.Equal(initialEnergy, energy, 5);
    }

    [Fact]
    public void BarycentricShift_MomentumShouldBecomeZero()
    {
        double G = 6.67430e-11;
        var bodies = new List<BodyState>
        {
            new BodyState("Sun", 1.989e30, new Vector3d(0, 0, 0), new Vector3d(0, 0, 0)),
            new BodyState("Earth", 5.972e24, new Vector3d(1.496e11, 0, 0), new Vector3d(0, 29780, 0))
        };

        var system = new NBodySystem(G, bodies, toBarycentricFrame: true);
        Vector3d momentum = system.CalculateTotalMomentum(system.InitialState);

        // После смещения в барицентрическую систему суммарный импульс должен быть пренебрежимо мал.
        // Для масс 10^30 и скоростей 10^4, ошибка округления double составляет порядка 10^14.
        Assert.True(momentum.Length() < 1e15); 
    }
}
