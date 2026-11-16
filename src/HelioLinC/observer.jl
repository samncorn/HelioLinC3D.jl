module Observer

using SPICE
using StaticArrays
using AstroTime

using Downloads: download



""" takes, very sepcifically, mpc observatory code info (longitude (deg), parallax constants)
"""
function observatory_position(l, c, s)
    return SVector{3}(
        c*cosd(l),
        c*sind(l),
        s
    )
end

"HAS TO BE JULIAN DAYS IN TDB"
JDTDB2ET(jd) = (jd - SPICE.j2000()) * 86400.0

function get_observatory_position_heliocentric(jd_tdb, obs_ref)
    et = JDTDB2ET(jd_tdb)
    pos, _ = SPICE.spkpos("EARTH", et, "J2000", "NONE", "SUN")
    M      = SPICE.pxform("ITRF93", "J2000", et)

    return pos + M*obs_ref
end

function furnish_spice_mpc(cache_path)
    kernels = [
        (true, "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp", "de440.bsp"),
        (true, "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_2026_260806_2126_predict.bpc", "earth_2026_260806_2126_predict.bpc"),
        (true, "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_620120_260806.bpc", "earth_620120_260806.bpc"),
        (true, "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/earth_latest_high_prec.bpc", "earth_latest_high_prec.bpc"),
        (true, "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls", "naif0012.tls"),
        (true, "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/pck00010.tpc", "pck00010.pck"),
        (false, "https://www.minorplanetcenter.net/iau/lists/ObsCodes.html", "ObsCodes.html"),
    ]

    # Download kernels (not sure if we actually need all of these)
    mkpath(cache_path)
    for (kern, f, targ) in kernels
        path = joinpath(cache_path, targ)
        if !isfile(path)
            @info "downloading $targ..."
            download(f, path)
        end
        if kern
            furnsh(path)
        end
    end
end

function load_observatories(cache_path)
    observatories = Dict{String, SVector{3, Float64}}()
    for l in readlines(joinpath(cache_path, "ObsCodes.html"))[3:end-1]
        try
            observatories[l[1:3]] = SVector{3}(
                parse(Float64, l[7:13]),
                parse(Float64, l[14:21]),
                parse(Float64, l[22:30]),
            )
        catch
            continue
        end
    end
    return observatories
end

function get_observer_position!(jd, cache)
    pos = if haskey(cache, jd)
        cache[jd]
    else
        _pos = SVector{3}(hl.Observer.get_observatory_position_heliocentric(jd, obs_ref)) / hl.AU2KM
        cache[jd] = _pos
        _pos
    end
    return pos
end

# end module
end