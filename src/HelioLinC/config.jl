#= structures for configuration =#
struct Hypothesis{T}
    gm  :: T
    r   :: T
    dr  :: T
    ddr :: T
    t   :: T
end

struct Arrow_Limits{T}
    r_min  :: T
    e_min  :: T
    e_max  :: T
    i_min  :: T
    i_max  :: T
    dt_min :: T
    dt_max :: T
end

# dbscan specific args
struct Cluster_Params{T}
    cluster_radius :: T
    max_rms        :: T
    vel_scale      :: T
    min_pts        :: Int
    max_pts        :: Int
    n_threads      :: Int
end
