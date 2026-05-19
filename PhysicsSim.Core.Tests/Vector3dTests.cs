using PhysicsSim.Core;
using Xunit;

namespace PhysicsSim.Core.Tests;

public class Vector3dTests
{
    [Fact]
    public void Zero_ReturnsAllComponents()
    {
        var v = Vector3d.Zero;
        Assert.Equal(0, v.X);
        Assert.Equal(0, v.Y);
        Assert.Equal(0, v.Z);
    }

    [Fact]
    public void Length_UnitVectors()
    {
        Assert.Equal(1.0, new Vector3d(1, 0, 0).Length(), 15);
        Assert.Equal(1.0, new Vector3d(0, 1, 0).Length(), 15);
        Assert.Equal(1.0, new Vector3d(0, 0, 1).Length(), 15);
    }

    [Fact]
    public void Length_ArbitraryVector()
    {
        var v = new Vector3d(3, 4, 0);
        Assert.Equal(5.0, v.Length(), 14);
    }

    [Fact]
    public void LengthSquared_MatchesLengthTimesLength()
    {
        var v = new Vector3d(2, 3, 6);
        Assert.Equal(v.Length() * v.Length(), v.LengthSquared(), 10);
    }

    [Fact]
    public void Normalized_UnitLength()
    {
        var v = new Vector3d(3, 4, 0).Normalized();
        Assert.Equal(1.0, v.Length(), 14);
        Assert.Equal(0.6, v.X, 14);
        Assert.Equal(0.8, v.Y, 14);
    }

    [Fact]
    public void Normalized_ZeroVector_ReturnsZero()
    {
        var v = Vector3d.Zero.Normalized();
        Assert.Equal(0, v.X);
        Assert.Equal(0, v.Y);
        Assert.Equal(0, v.Z);
    }

    [Fact]
    public void Dot_OrthogonalVectors_IsZero()
    {
        var a = new Vector3d(1, 0, 0);
        var b = new Vector3d(0, 1, 0);
        Assert.Equal(0, Vector3d.Dot(a, b), 15);
    }

    [Fact]
    public void Dot_ParallelVectors()
    {
        var a = new Vector3d(2, 3, 4);
        var b = new Vector3d(2, 3, 4);
        Assert.Equal(a.LengthSquared(), Vector3d.Dot(a, b), 10);
    }

    [Fact]
    public void Cross_UnitVectors_RightHandRule()
    {
        var x = new Vector3d(1, 0, 0);
        var y = new Vector3d(0, 1, 0);
        var z = Vector3d.Cross(x, y);
        Assert.Equal(0, z.X, 15);
        Assert.Equal(0, z.Y, 15);
        Assert.Equal(1, z.Z, 15);
    }

    [Fact]
    public void Cross_ParallelVectors_IsZero()
    {
        var a = new Vector3d(1, 2, 3);
        var b = new Vector3d(2, 4, 6);
        var c = Vector3d.Cross(a, b);
        Assert.True(c.Length() < 1e-14);
    }

    [Fact]
    public void Addition()
    {
        var a = new Vector3d(1, 2, 3);
        var b = new Vector3d(4, 5, 6);
        var c = a + b;
        Assert.Equal(5, c.X);
        Assert.Equal(7, c.Y);
        Assert.Equal(9, c.Z);
    }

    [Fact]
    public void Subtraction()
    {
        var a = new Vector3d(5, 7, 9);
        var b = new Vector3d(1, 2, 3);
        var c = a - b;
        Assert.Equal(4, c.X);
        Assert.Equal(5, c.Y);
        Assert.Equal(6, c.Z);
    }

    [Fact]
    public void ScalarMultiplication_LeftAndRight()
    {
        var v = new Vector3d(1, 2, 3);
        var left = 2.0 * v;
        var right = v * 2.0;
        Assert.Equal(left.X, right.X);
        Assert.Equal(left.Y, right.Y);
        Assert.Equal(left.Z, right.Z);
        Assert.Equal(2, left.X);
        Assert.Equal(4, left.Y);
        Assert.Equal(6, left.Z);
    }

    [Fact]
    public void Division()
    {
        var v = new Vector3d(6, 8, 10);
        var r = v / 2.0;
        Assert.Equal(3, r.X);
        Assert.Equal(4, r.Y);
        Assert.Equal(5, r.Z);
    }

    [Fact]
    public void Negation()
    {
        var v = new Vector3d(1, -2, 3);
        var n = -v;
        Assert.Equal(-1, n.X);
        Assert.Equal(2, n.Y);
        Assert.Equal(-3, n.Z);
    }

    [Fact]
    public void RecordEquality()
    {
        var a = new Vector3d(1, 2, 3);
        var b = new Vector3d(1, 2, 3);
        Assert.Equal(a, b);
    }
}
