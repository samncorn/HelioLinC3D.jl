### Installation
- Install the [Julia runtime](https://julialang.org/downloads/)
- Download this package
    ``` 
    git clone https://github.com/samncorn/HelioLinC3D.jl.git
    ```

## Basic Usage

### Clustering

The `heliolinc.jl` program ingests optical measurments and a set of search parameters, generating a set of files with candidate linkages.

```
julia -t [n_threads] path_to_project/scripts/heliolinc.jl [args...]
```
where path_to_project is the path to the downloaded repository.

When running for the first time, the clustering program will download ancillary data to a .cache directory in the project directory automatically.

### Orbit filtering

The candidate linkages typically have a large amount of overlap and junk, requiring further processing with the `filter-linkages.jl` program, run in the same manner as the `heliolinc.jl` program. 

```
julia -t [n_threads] path_to_project/scripts/filter-linkages.jl [args...]
```

## Basic Usage

### Clustering

To Do

### Orbit filtering

To Do

Many of the named arguments are different, but if using a config file the same one as
for `heliolinc.jl` will work, and is recommended.