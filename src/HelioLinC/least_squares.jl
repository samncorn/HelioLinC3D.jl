"""
f(x) is a function taking x that computes the residuals for each measurement
f(x) and df(x) must be iterators which match the measurements length

so if each measurement is a scalar, each iterate of df is the gradient of f 
if each measurement is a vector, each iterate of df is the jacobian

if the first case, needs to a 1 row matrix, or a transposed vector

NOTE : f(x) and df(x) can be given as lazy iterators to minimize allocations

TODO: improve the gradient descent coefficient determination (add an argument which takes a function of the step, returning the change in gradient descent  coefficient)
"""
function _levenberg_marquardt(f, df, x0; mu_0 = 100.0, mu_f = 5.0, rel_tol=1e-10, abs_tol = 1e-10, max_iter=100, sink = nothing)
    x  = x0
    mu = mu_0
    I  = Diagonal(zero(x) .+ 1.0)
    # initial cost
    dy = f(x)
    Ji = resid_cost(dy)
    i  = 0
    rel_err = Inf
    abs_err = Inf
    while rel_err > rel_tol && abs_err > abs_tol && i < max_iter
        dy = f(x)
        H  = df(x)
        # if appropriate inputs are passed, these should end up being sized matrices/vectors
        HtH  = accumulate_HtH(H)
        Htdy = accumulate_Htdy(H, dy)
        dx   = pinv(HtH + mu*I)*Htdy
        # update residuals
        dy = f(x - dx)
        J  = resid_cost(dy)

        # check how well our step went. If bad, increasing gradient descent component
        # and redo step
        if J > Ji
            mu *= mu_f
            i += 1
            continue
        end

        if !isnothing(sink)
            push!(sink, (i = i, x = x, dx = -dx, J = J, mu = mu))
        end

        mu /= mu_f
        x -= dx
        rel_err = abs((Ji - J)/Ji)
        abs_err = abs(Ji - J)
        Ji  = J
        i  += 1
    end
    if !isnothing(sink)
        push!(sink, (i = i, x = x, dx = 0.0, J = Ji, mu = mu))
    end
    return x, dy, i
end

function levenberg_marquardt(f, df, x0; mu_0 = 100.0, mu_f = 5.0, rel_tol=1e-10, abs_tol = 1e-10, max_iter = 100, sink = nothing)
    x, dy, i = _levenberg_marquardt(f, df, x0; mu_0 = 100.0, mu_f = 5.0, rel_tol = rel_tol, abs_tol = abs_tol, max_iter = max_iter, sink = sink)
    i >= max_iter && throw(MaximumIterationError(i))
    return x, dy
end

# convenient functions which will exploit iterators
resid_cost(resid)      = reduce(+, Iterators.map(x -> x^2, Iterators.flatten(resid)))/2
accumulate_HtH(H)      = reduce(+, Iterators.map(h -> transpose(h)*h, H))
accumulate_Htdy(H, dy) = reduce(+, Iterators.map(x -> transpose(x[1])*x[2], zip(H, dy)))
