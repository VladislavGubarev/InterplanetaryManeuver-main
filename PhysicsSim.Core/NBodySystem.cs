using System.Numerics;

namespace PhysicsSim.Core;

/// <summary>
/// Ньютоновская N-body модель.
/// Для каждого тела хранятся координаты и скорости:
/// r_i' = v_i
/// v_i' = sum_{j!=i} G * m_j * (r_j - r_i) / |r_j - r_i|^3
/// </summary>
public sealed class NBodySystem : IOdeSystem
{
    public double GravitationalConstant { get; }
    public string[] Names { get; }
    public double[] Masses { get; }

    public int BodyCount => Masses.Length;
    public int Dimension => BodyCount * 6; // На тело: x, y, z, vx, vy, vz.

    public NBodySystem(double gravitationalConstant, IReadOnlyList<BodyState> bodies, bool toBarycentricFrame = true)
    {
        if (bodies.Count < 2)
            throw new ArgumentException("Need at least 2 bodies.", nameof(bodies));

        GravitationalConstant = gravitationalConstant;
        Names = new string[bodies.Count];
        Masses = new double[bodies.Count];

        for (int i = 0; i < bodies.Count; i++)
        {
            Names[i] = bodies[i].Name;
            Masses[i] = bodies[i].Mass;
        }

        InitialState = new double[Dimension];
        for (int i = 0; i < bodies.Count; i++)
            WriteBodyState(InitialState, i, bodies[i].Position, bodies[i].Velocity);

        if (toBarycentricFrame)
            ShiftInitialStateToBarycentricFrame();
    }

    /// <summary>
    /// Начальный вектор состояния.
    /// Формат: [x, y, z, vx, vy, vz] для каждого тела подряд.
    /// </summary>
    public double[] InitialState { get; }

    public void ComputeDerivatives(double t, ReadOnlySpan<double> y, Span<double> dy)
    {
        if (y.Length != Dimension) throw new ArgumentException("Invalid state length.", nameof(y));
        if (dy.Length != Dimension) throw new ArgumentException("Invalid derivative length.", nameof(dy));

        // Кинематика: производная координаты равна текущей скорости.
        for (int i = 0; i < BodyCount; i++)
        {
            int baseIdx = i * 6;
            dy[baseIdx + 0] = y[baseIdx + 3];
            dy[baseIdx + 1] = y[baseIdx + 4];
            dy[baseIdx + 2] = y[baseIdx + 5];
        }

        // Динамика: суммируем гравитационные ускорения.
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

                double mj = Masses[j];
                if (mj == 0) continue;

                int bj = j * 6;
                double dx = y[bj + 0] - xi;
                double dyv = y[bj + 1] - yi;
                double dz = y[bj + 2] - zi;

                double r2 = dx * dx + dyv * dyv + dz * dz;
                if (r2 <= 1e-18) continue;

                // Ускорение a = G * mj / r^3 * r_vec
                double invR3 = Math.Pow(r2, -1.5);
                double k = GravitationalConstant * mj * invR3;

                ax += k * dx;
                ay += k * dyv;
                az += k * dz;
            }

            dy[bi + 3] = ax;
            dy[bi + 4] = ay;
            dy[bi + 5] = az;
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
                double mj = Masses[j];
                double dx = state[bj + 0] - state[bi + 0];
                double dy = state[bj + 1] - state[bi + 1];
                double dz = state[bj + 2] - state[bi + 2];
                double r2 = dx * dx + dy * dy + dz * dz;

                if (r2 > 1e-18)
                {
                    totalPotential -= GravitationalConstant * mi * mj / Math.Sqrt(r2);
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
