function observer_interpolant(t, x, y, z, dx, dy, dz)
    x_interp = CubicHermiteSpline(dx, x, t)
    y_interp = CubicHermiteSpline(dy, y, t)
    z_interp = CubicHermiteSpline(dz, z, t)
    return s -> SVector{3}((x_interp(s), y_interp(s), z_interp(s)))
end

function cubic_hermite_spline(t, x, y, z, dx, dy, dz)
    x_interp = CubicHermiteSpline(dx, x, t)
    y_interp = CubicHermiteSpline(dy, y, t)
    z_interp = CubicHermiteSpline(dz, z, t)
    return s -> SVector{3}((x_interp(s), y_interp(s), z_interp(s)))
end

function cubic_spline(t, x, y, z)
    x_interp = CubicSpline(x, t)
    y_interp = CubicSpline(y, t)
    z_interp = CubicSpline(z, t)
    return s -> SVector{3}((x_interp(s), y_interp(s), z_interp(s)))
end