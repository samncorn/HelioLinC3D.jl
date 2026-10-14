# solar_gravity(pos, gm) = -gm/norm(pos)^3 * pos

""" Solve the lambert problem numerically (shooting method) for the single revolution solution
"""
function _lambert_dns(pos1, pos2, dt, propagate; max_iter = 20, abs_tol = 1e-12, rel_tol = 1e-12)
    # finite diff the velocity for initial guess
    vel1 = (pos2 - pos1)/dt
    iter = 1
    rel_err  = Inf
    abs_err  = Inf
    while iter < max_iter && abs_err > abs_tol && rel_err > rel_tol
        posf, _ = propagate(pos1, vel1, dt)
        dvel  = (pos2 - posf)/dt
        vel1 += dvel
        abs_err = norm(dvel)
        rel_err = norm(dvel)/norm(vel1)
        iter += 1
    end
    return iter, vel1
end

function lambert_dns(pos1, pos2, dt, propagate; max_iter = 20, abs_tol = 1e-12, rel_tol = 1e-12)
    i, vel = _lambert_dns(pos1, pos2, dt, propagate; max_iter = max_iter, abs_tol = abs_tol, rel_tol = rel_tol)
    i >= max_iter && throw(MaximumIterationError(i))
    return vel
end

"""
terminates after either meeting tolerance or reacing max_iter. Throws an error for max_iter
"""
function _apparent_state(target_traj, observer, t, c; lt_tol = 1e-15, max_iter = 20)
    # initialize
    app_state = target_traj(t) - observer
    lti = norm(app_state) / c
    @debug("$(app_state)")
    dlt = Inf
    i   = 0
    while dlt > lt_tol && i < max_iter
        app_state = target_traj(t - lti) - observer
        lt  = norm(app_state) / c
        dlt = abs(lt - lti)
        @debug("$(app_state), $(lt), $(dlt)")
        lti = lt
        i += 1
    end
    return target_traj(t - lti) - observer, t-lti, i
end

function apparent_state(target_traj, observer, t, c; lt_tol = 1e-12, max_iter = 20)
    state, lt, i = _apparent_state(target_traj, observer, t, c; lt_tol = lt_tol, max_iter = max_iter)
    if i >= max_iter
        @debug("apparent state failed after $(i) iters with tol $(lt_tol)")
        @debug("t: $(t), dlt: $(t - lt)")
        @debug("c: $(c)")
        throw(MaximumIterationError(i))
    end
    return state
end

function angles_residual(ang1, ang2)
    dang = ang2 - ang1
    # basis vectors
    v_xy = normalize(SVector{3}((ang1[1], ang1[2], 0.0)))
    zhat = SVector{3}((0.0, 0.0, 1.0))
    A = normalize(-cross(v_xy, zhat))
    D = normalize(cross(ang1, A))
    return dot(dang, A), dot(dang, D)
end

"return an iterator of residuals"
function orbit_residuals(measurements, traj, c)
    resids = Iterators.map(
        m -> SVector{2}(angles_residual(
            m.ang, 
            normalize(_apparent_state(
                traj, 
                m.pos, 
                m.time, c;
                max_iter = 3
                )[1])
            )), 
        measurements
        )
    return resids
end

""" generate an interpolant from t1 to t2 given initial position and dynamics
dynamics should be a function on the pos-vel state
requires t1 and t2 be ordered properly
"""
function integrate_traj_interpolant(pos0, vel0, t1, t2, dynamics; tol = 1e-15)
    @assert t1 < t2
    ode = ODEProblem{false}(dynamics, SA[pos0..., vel0...], (t1, t2))
    sol = solve(ode, Tsit5(), reltol = tol, abstol = tol)
    return sol
end

""" Requires an interpolant constructing function to properly handle light times. This is uncoupled from the propagate
function used to solve the short arc lambert problem, but the user should specify an interpant contructor that uses the
given propagate.
"""
function herget_residuals(measurements, mi, mj, rhoi, rhoj, c, propagate, interpolant_constructor)
    ti   = mi.time - rhoi/c
    tj   = mj.time - rhoj/c
    posi = mi.pos + mi.ang*rhoi
    posj = mj.pos + mj.ang*rhoj
    veli = lambert_dns(posi, posj, tj - ti, propagate; max_iter = 20, abs_tol = 1e-12, rel_tol = 1e-12)
    # now we have an estimate of the orbit. Build an interpolator across the trajectory
    times    = Iterators.map(m -> m.time, measurements)
    min_time = minimum(times)
    max_time = maximum(times)
    traj     = interpolant_constructor(posi, veli, min_time, max_time) # propagate(posi, veli, times)
    return orbit_residuals(measurements, traj, c)
end

pq_rho_partials(PQ11, PQ12, PQ21, PQ22, drho) = SMatrix{2, 2}((
    PQ11[1] - PQ12[1], # p
    PQ11[2] - PQ12[2], # q
    PQ21[1] - PQ22[1], # p
    PQ21[2] - PQ22[2], # q 
    )) ./ (2drho)

""" central difference the hell out of it
"""
function herget_partials(measurements, mi, mj, rhoi, rhoj, c, propagate, interpolant_constructor; drho = 0.001)
    PQ11 = herget_residuals(measurements, mi, mj, rhoi + drho, rhoj,        c, propagate, interpolant_constructor)
    PQ12 = herget_residuals(measurements, mi, mj, rhoi - drho, rhoj,        c, propagate, interpolant_constructor)
    PQ21 = herget_residuals(measurements, mi, mj, rhoi       , rhoj + drho, c, propagate, interpolant_constructor)
    PQ22 = herget_residuals(measurements, mi, mj, rhoi       , rhoj - drho, c, propagate, interpolant_constructor)
    return map(pq -> pq_rho_partials(pq..., drho), zip(PQ11, PQ12, PQ21, PQ22))
end

function herget_od(measurements, i, j, rho0i, rho0j, c, propagate, interpolant_constructor; max_iter = 20, rel_tol = 1e-12, abs_tol = 1e-8, drho = 0.001, log_sink = nothing, mu_0 = 100.0)
    # levenberg-marquardt the hell out of this
    # takes the residual function, not the measurement function
    idx = collect(Iterators.filter(x -> x != i && x != j, eachindex(measurements)))
    m_subset = @view measurements[idx]

    x   = SVector{2}((rho0i, rho0j))
    f  = x -> herget_residuals(m_subset, measurements[i], measurements[j], x..., c, propagate, interpolant_constructor) 
    df = x -> herget_partials(m_subset, measurements[i], measurements[j], x..., c, propagate, interpolant_constructor; drho = drho)
    (rhoi, rhoj), resids, n_iters = _levenberg_marquardt(f, df, x; sink = log_sink, rel_tol = rel_tol, abs_tol = abs_tol, mu_0 = mu_0);

    converged = n_iters <= max_iter
    # the determiend state
    ti = measurements[i].time
    tj = measurements[j].time
    posi = measurements[i].pos + rhoi*measurements[i].ang
    posj = measurements[j].pos + rhoj*measurements[j].ang
    veli = lambert_dns(posi, posj, tj - ti, propagate)
    state = SA[posi..., veli...]
    # return OD_Solution(converged, state, resid, ti)
    return converged, n_iters, state, ti, rhoi, rhoj, resids
end