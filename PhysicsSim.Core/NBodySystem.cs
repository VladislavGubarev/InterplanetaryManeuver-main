using System.Numerics;
using System.Numerics.Tensors;

namespace PhysicsSim.Core;

/// <summary>
/// Ньютоновская N-body модель.
/// Для каждого тела хранятся координаты и скорости:
/// r_i' = v_i
/// v_i' = sum_{j!=i} GM_j * (r_j - r_i) / |r_j - r_i|^3
///
/// Используем гравитационные параметры GM вместо G*M — это стандарт
/// в аэрокосмической отрасли, т.к. GM известен с точностью ~1e-10,
/// тогда как G и M по отдельности — ~1e-5.
/// </summary>
public sealed class NBodySystem : IOdeSystem
{
    /// <summary>Softening по умолчанию: ε² предотвращает деление на ноль
    /// при очень малых расстояниях без потери точности на планетарных масштабах.
    /// ε = 100 м → ε² = 10⁴ м² — пренебрежимо мало на расстояниях > 10⁶ м.</summary>
    public const double DefaultSofteningSquared = 1e4;

    private readonly double _softeningSquared;
    private readonly bool _useTensors;

    // Предвыделенные буферы для тензорного пути (SoA-layout)
    private readonly double[]? _px, _py, _pz;
    private readonly double[]? _diffX, _diffY, _diffZ;
    private readonly double[]? _r2Buf, _invR3Buf, _forceBuf;

    public string[] Names { get; }
    public double[] Masses { get; }
    public double[] GMs { get; }
    public bool UseTensors => _useTensors;

    public int BodyCount => GMs.Length;
    public int Dimension => BodyCount * 6;

    public NBodySystem(double gravitationalConstant, IReadOnlyList<BodyState> bodies,
        bool toBarycentricFrame = true, double softeningSquared = DefaultSofteningSquared,
        bool useTensors = false)
    {
        if (bodies.Count < 2)
            throw new ArgumentException("Need at least 2 bodies.", nameof(bodies));

        _softeningSquared = softeningSquared;
        _useTensors = useTensors;
        Names = new string[bodies.Count];
        Masses = new double[bodies.Count];
        GMs = new double[bodies.Count];

        for (int i = 0; i < bodies.Count; i++)
        {
            Names[i] = bodies[i].Name;
            Masses[i] = bodies[i].Mass;
            GMs[i] = gravitationalConstant * bodies[i].Mass;
        }

        InitialState = new double[Dimension];
        for (int i = 0; i < bodies.Count; i++)
            WriteBodyState(InitialState, i, bodies[i].Position, bodies[i].Velocity);

        if (toBarycentricFrame)
            ShiftInitialStateToBarycentricFrame();

        if (_useTensors)
            AllocateTensorBuffers(bodies.Count, out _px, out _py, out _pz, out _diffX, out _diffY, out _diffZ, out _r2Buf, out _invR3Buf, out _forceBuf);
    }

    public NBodySystem(IReadOnlyList<BodyState> bodies, IReadOnlyList<double> gms,
        bool toBarycentricFrame = true, double softeningSquared = DefaultSofteningSquared,
        bool useTensors = false)
    {
        if (bodies.Count < 2)
            throw new ArgumentException("Need at least 2 bodies.", nameof(bodies));
        if (bodies.Count != gms.Count)
            throw new ArgumentException("GM count must match body count.", nameof(gms));

        _softeningSquared = softeningSquared;
        _useTensors = useTensors;
        Names = new string[bodies.Count];
        Masses = new double[bodies.Count];
        GMs = new double[bodies.Count];

        for (int i = 0; i < bodies.Count; i++)
        {
            Names[i] = bodies[i].Name;
            Masses[i] = bodies[i].Mass;
            GMs[i] = gms[i];
        }

        InitialState = new double[Dimension];
        for (int i = 0; i < bodies.Count; i++)
            WriteBodyState(InitialState, i, bodies[i].Position, bodies[i].Velocity);

        if (toBarycentricFrame)
            ShiftInitialStateToBarycentricFrame();

        if (_useTensors)
            AllocateTensorBuffers(bodies.Count, out _px, out _py, out _pz, out _diffX, out _diffY, out _diffZ, out _r2Buf, out _invR3Buf, out _forceBuf);
    }

    private static void AllocateTensorBuffers(int n,
        out double[] px, out double[] py, out double[] pz,
        out double[] diffX, out double[] diffY, out double[] diffZ,
        out double[] r2Buf, out double[] invR3Buf, out double[] forceBuf)
    {
        px = new double[n];
        py = new double[n];
        pz = new double[n];
        diffX = new double[n];
        diffY = new double[n];
        diffZ = new double[n];
        r2Buf = new double[n];
        invR3Buf = new double[n];
        forceBuf = new double[n];
    }

    public double[] InitialState { get; }

    public void ComputeDerivatives(double t, ReadOnlySpan<double> y, Span<double> dy)
    {
        if (y.Length != Dimension) throw new ArgumentException("Invalid state length.", nameof(y));
        if (dy.Length != Dimension) throw new ArgumentException("Invalid derivative length.", nameof(dy));

        for (int i = 0; i < BodyCount; i++)
        {
            int baseIdx = i * 6;
            dy[baseIdx + 0] = y[baseIdx + 3];
            dy[baseIdx + 1] = y[baseIdx + 4];
            dy[baseIdx + 2] = y[baseIdx + 5];
        }

        if (_useTensors)
            ComputeAccelerationsTensor(y, dy);
        else
            ComputeAccelerationsScalar(y, dy);
    }

    private void ComputeAccelerationsScalar(ReadOnlySpan<double> y, Span<double> dy)
    {
        for (int i = 0; i < BodyCount; i++)
        {
            int bi = i * 6;
            double xi = y[bi + 0];
            double yi = y[bi + 1];
            double zi = y[bi + 2];

            double ax = 0, ay = 0, az = 0;

            for (int j = 0; j < BodyCount; j++)
            {
                if (j == i) continue;

                double gmj = GMs[j];
                if (gmj == 0) continue;

                int bj = j * 6;
                double dx = y[bj + 0] - xi;
                double dyv = y[bj + 1] - yi;
                double dz = y[bj + 2] - zi;

                double r2 = dx * dx + dyv * dyv + dz * dz + _softeningSquared;

                double invR3 = Math.Pow(r2, -1.5);
                double k = gmj * invR3;

                ax += k * dx;
                ay += k * dyv;
                az += k * dz;
            }

            dy[bi + 3] = ax;
            dy[bi + 4] = ay;
            dy[bi + 5] = az;
        }
    }

    /// <summary>
    /// Тензорный путь: использует TensorPrimitives (SIMD) для вычисления
    /// гравитационных ускорений. Для каждого тела i все разности координат,
    /// расстояния и силы вычисляются параллельно по всем j-телам через
    /// векторизованные операции над Span&lt;double&gt;.
    /// </summary>
    private void ComputeAccelerationsTensor(ReadOnlySpan<double> y, Span<double> dy)
    {
        int n = BodyCount;

        // Распаковка позиций из AoS (interleaved) в SoA (раздельные массивы)
        for (int j = 0; j < n; j++)
        {
            int bj = j * 6;
            _px![j] = y[bj + 0];
            _py![j] = y[bj + 1];
            _pz![j] = y[bj + 2];
        }

        Span<double> diffX = _diffX.AsSpan(0, n);
        Span<double> diffY = _diffY.AsSpan(0, n);
        Span<double> diffZ = _diffZ.AsSpan(0, n);
        Span<double> r2 = _r2Buf.AsSpan(0, n);
        Span<double> invR3 = _invR3Buf.AsSpan(0, n);
        Span<double> force = _forceBuf.AsSpan(0, n);

        ReadOnlySpan<double> pxSpan = _px.AsSpan(0, n);
        ReadOnlySpan<double> pySpan = _py.AsSpan(0, n);
        ReadOnlySpan<double> pzSpan = _pz.AsSpan(0, n);
        ReadOnlySpan<double> gmSpan = GMs.AsSpan(0, n);

        for (int i = 0; i < n; i++)
        {
            int bi = i * 6;
            double xi = y[bi + 0];
            double yi_val = y[bi + 1];
            double zi = y[bi + 2];

            // diffX[j] = px[j] - xi  (для всех j одновременно через SIMD)
            TensorPrimitives.Subtract(pxSpan, xi, diffX);
            TensorPrimitives.Subtract(pySpan, yi_val, diffY);
            TensorPrimitives.Subtract(pzSpan, zi, diffZ);

            // r2[j] = diffX[j]² + diffY[j]² + diffZ[j]² + ε²
            TensorPrimitives.Multiply(diffX, diffX, r2);
            TensorPrimitives.MultiplyAdd((ReadOnlySpan<double>)diffY, (ReadOnlySpan<double>)diffY, r2, r2);
            TensorPrimitives.MultiplyAdd((ReadOnlySpan<double>)diffZ, (ReadOnlySpan<double>)diffZ, r2, r2);
            TensorPrimitives.Add(r2, _softeningSquared, r2);

            // invR3[j] = r2[j]^(-1.5)
            for (int j = 0; j < n; j++)
                invR3[j] = Math.Pow(r2[j], -1.5);

            // force[j] = GM[j] * invR3[j]
            TensorPrimitives.Multiply(gmSpan, (ReadOnlySpan<double>)invR3, force);

            // Обнуляем вклад i-го тела на себя
            force[i] = 0;

            // ax = Σ force[j] * diffX[j],  ay = ...,  az = ...
            TensorPrimitives.Multiply((ReadOnlySpan<double>)force, (ReadOnlySpan<double>)diffX, diffX);
            TensorPrimitives.Multiply((ReadOnlySpan<double>)force, (ReadOnlySpan<double>)diffY, diffY);
            TensorPrimitives.Multiply((ReadOnlySpan<double>)force, (ReadOnlySpan<double>)diffZ, diffZ);

            dy[bi + 3] = TensorPrimitives.Sum((ReadOnlySpan<double>)diffX);
            dy[bi + 4] = TensorPrimitives.Sum((ReadOnlySpan<double>)diffY);
            dy[bi + 5] = TensorPrimitives.Sum((ReadOnlySpan<double>)diffZ);
        }
    }

    public BodyState ReadBodyFromState(ReadOnlySpan<double> state, int bodyIndex)
    {
        int b = bodyIndex * 6;
        var position = new Vector3d(state[b + 0], state[b + 1], state[b + 2]);
        var velocity = new Vector3d(state[b + 3], state[b + 4], state[b + 5]);
        return new BodyState(Names[bodyIndex], Masses[bodyIndex], position, velocity);
    }

    public static void WriteBodyState(Span<double> state, int bodyIndex, Vector3d position, Vector3d velocity)
    {
        int b = bodyIndex * 6;
        state[b + 0] = position.X;
        state[b + 1] = position.Y;
        state[b + 2] = position.Z;
        state[b + 3] = velocity.X;
        state[b + 4] = velocity.Y;
        state[b + 5] = velocity.Z;
    }

    public void ShiftInitialStateToBarycentricFrame()
    {
        var comPos = Vector3d.Zero;
        var comVel = Vector3d.Zero;
        double mTotal = 0;

        for (int i = 0; i < BodyCount; i++)
        {
            double m = Masses[i];
            if (m <= 0) continue;

            mTotal += m;
            int b = i * 6;
            comPos += m * new Vector3d(InitialState[b + 0], InitialState[b + 1], InitialState[b + 2]);
            comVel += m * new Vector3d(InitialState[b + 3], InitialState[b + 4], InitialState[b + 5]);
        }

        if (mTotal == 0) return;

        comPos /= mTotal;
        comVel /= mTotal;

        for (int i = 0; i < BodyCount; i++)
        {
            int b = i * 6;
            InitialState[b + 0] -= comPos.X;
            InitialState[b + 1] -= comPos.Y;
            InitialState[b + 2] -= comPos.Z;
            InitialState[b + 3] -= comVel.X;
            InitialState[b + 4] -= comVel.Y;
            InitialState[b + 5] -= comVel.Z;
        }
    }

    public double CalculateTotalEnergy(ReadOnlySpan<double> state)
    {
        double totalKinetic = 0;
        double totalPotential = 0;

        for (int i = 0; i < BodyCount; i++)
        {
            int bi = i * 6;
            double mi = Masses[i];
            double vx = state[bi + 3];
            double vy = state[bi + 4];
            double vz = state[bi + 5];
            double v2 = vx * vx + vy * vy + vz * vz;
            totalKinetic += 0.5 * mi * v2;

            for (int j = i + 1; j < BodyCount; j++)
            {
                int bj = j * 6;
                double dx = state[bj + 0] - state[bi + 0];
                double dy = state[bj + 1] - state[bi + 1];
                double dz = state[bj + 2] - state[bi + 2];
                double r2 = dx * dx + dy * dy + dz * dz;

                if (r2 > 1e-18)
                {
                    // Для энергии используем точный GM вместо G*mi*mj:
                    // E_pot = -GM_i * m_j / r = -GM_j * m_i / r (эквивалентно).
                    totalPotential -= GMs[i] * Masses[j] / Math.Sqrt(r2);
                }
            }
        }
        return totalKinetic + totalPotential;
    }

    public Vector3d CalculateTotalMomentum(ReadOnlySpan<double> state)
    {
        double px = 0, py = 0, pz = 0;
        for (int i = 0; i < BodyCount; i++)
        {
            int bi = i * 6;
            double mi = Masses[i];
            px += mi * state[bi + 3];
            py += mi * state[bi + 4];
            pz += mi * state[bi + 5];
        }
        return new Vector3d(px, py, pz);
    }
}
