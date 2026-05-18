using PhysicsSim.Core;
using PhysicsSim.Core.Ode;
using Xunit;

namespace PhysicsSim.Core.Tests;

public class DormandPrince45Tests
{
    /// <summary>
    /// Simple harmonic oscillator: y'' = -y
    /// Encoded as y[0]=position, y[1]=velocity.
    /// </summary>
    private sealed class HarmonicOscillator : IOdeSystem
    {
        public int Dimension => 2;

        public void ComputeDerivatives(double t, ReadOnlySpan<double> y, Span<double> dy)
        {
            dy[0] = y[1];
            dy[1] = -y[0];
        }
    }

    /// <summary>
    /// Exponential decay: y' = -y
    /// Analytic solution: y(t) = y0 * exp(-t)
    /// </summary>
    private sealed class ExponentialDecay : IOdeSystem
    {
        public int Dimension => 1;

        public void ComputeDerivatives(double t, ReadOnlySpan<double> y, Span<double> dy)
        {
            dy[0] = -y[0];
        }
    }

    [Fact]
    public void HarmonicOscillator_PreservesAmplitude()
    {
        var system = new HarmonicOscillator();
        var settings = new IntegrationSettings
        {
            AbsTol = 1e-10,
            RelTol = 1e-10,
            InitialStep = 0.01,
            MinStep = 1e-12,
            MaxStep = 1.0,
        };

        double[] y0 = [1.0, 0.0]; // cos(0)=1, sin(0)=0
        double t0 = 0, t1 = 2 * Math.PI;
        double outputDt = 0.1;

        var samples = new List<(double t, double y, double v)>();
        DormandPrince45.IntegrateFixedOutput(
            system, t0, y0, t1, outputDt, settings,
            (t, y) => samples.Add((t, y[0], y[1])));

        Assert.True(samples.Count > 10);

        var last = samples[^1];
        Assert.Equal(1.0, last.y, 6);
        Assert.Equal(0.0, last.v, 5);
    }

    [Fact]
    public void ExponentialDecay_MatchesAnalyticSolution()
    {
        var system = new ExponentialDecay();
        var settings = new IntegrationSettings
        {
            AbsTol = 1e-12,
            RelTol = 1e-12,
            InitialStep = 0.01,
            MinStep = 1e-15,
            MaxStep = 1.0,
        };

        double y0Val = 5.0;
        double[] y0 = [y0Val];
        double t0 = 0, t1 = 3.0;
        double outputDt = 0.5;

        var samples = new List<(double t, double y)>();
        DormandPrince45.IntegrateFixedOutput(
            system, t0, y0, t1, outputDt, settings,
            (t, y) => samples.Add((t, y[0])));

        foreach (var (t, y) in samples)
        {
            double expected = y0Val * Math.Exp(-t);
            Assert.Equal(expected, y, 8);
        }
    }

    [Fact]
    public void EmitsFirstAndLastSample()
    {
        var system = new ExponentialDecay();
        var settings = new IntegrationSettings
        {
            AbsTol = 1e-8,
            RelTol = 1e-8,
            InitialStep = 0.01,
        };

        double[] y0 = [1.0];
        var times = new List<double>();

        DormandPrince45.IntegrateFixedOutput(
            system, 0, y0, 2.0, 0.5, settings,
            (t, _) => times.Add(t));

        Assert.Equal(0.0, times[0], 12);
        Assert.Equal(2.0, times[^1], 12);
    }

    [Fact]
    public void StopCondition_TerminatesEarly()
    {
        var system = new ExponentialDecay();
        var settings = new IntegrationSettings
        {
            AbsTol = 1e-8,
            RelTol = 1e-8,
            InitialStep = 0.01,
        };

        double[] y0 = [10.0];
        var times = new List<double>();

        string? reason = DormandPrince45.IntegrateFixedOutput(
            system, 0, y0, 100.0, 0.5, settings,
            (t, _) => times.Add(t),
            (t, y) => y[0] < 1.0 ? "Decayed below threshold" : null);

        Assert.NotNull(reason);
        Assert.Contains("threshold", reason);
        Assert.True(times[^1] < 100.0);
    }

    [Fact]
    public void CancellationToken_StopsIntegration()
    {
        var system = new ExponentialDecay();
        var settings = new IntegrationSettings
        {
            AbsTol = 1e-8,
            RelTol = 1e-8,
            InitialStep = 0.01,
        };

        var cts = new CancellationTokenSource();
        double[] y0 = [1.0];

        cts.Cancel();

        Assert.Throws<OperationCanceledException>(() =>
            DormandPrince45.IntegrateFixedOutput(
                system, 0, y0, 100.0, 0.1, settings,
                (_, _) => { },
                cancellationToken: cts.Token));
    }

    [Fact]
    public void InvalidArguments_Throw()
    {
        var system = new ExponentialDecay();
        var settings = new IntegrationSettings();
        double[] y0 = [1.0];

        Assert.Throws<ArgumentException>(() =>
            DormandPrince45.IntegrateFixedOutput(system, 5.0, y0, 3.0, 0.1, settings, (_, _) => { }));

        Assert.Throws<ArgumentException>(() =>
            DormandPrince45.IntegrateFixedOutput(system, 0, y0, 5.0, -0.1, settings, (_, _) => { }));

        double[] wrongDim = [1.0, 2.0];
        Assert.Throws<ArgumentException>(() =>
            DormandPrince45.IntegrateFixedOutput(system, 0, wrongDim, 5.0, 0.1, settings, (_, _) => { }));
    }
}
