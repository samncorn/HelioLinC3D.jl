############################################
# CONSTANTS
###########################################
const C_au_day = 173.145 # speed of light in au/day
const C_km_s = 299_792.458 # km/s

# Gaussian Gravitational Constant [au^1.5/Msun^0.5/D]
const GAUSS_K = 0.01720209894846

# default gravitational parameter [Gaussian units]
const GM_sun = GAUSS_K*GAUSS_K

# Julian date of J2000 epoch
const J2000_JD = 2451545.0

# Earth obliquity [deg]
const EARTH_OBLIQUITY = 23.4392911111112

const AU2KM = 149597870.7

I3 = SMatrix{3, 3, Float64}([
    1.0 0.0 0.0;
    0.0 1.0 0.0;
    0.0 0.0 1.0;
])

# ecliptic to icrf and back
const ICRF2ECL = SMatrix{3, 3, Float64}([
    1.0  0.0                   0.0                  ;
    0.0  cosd(EARTH_OBLIQUITY) sind(EARTH_OBLIQUITY);
    0.0 -sind(EARTH_OBLIQUITY) cosd(EARTH_OBLIQUITY);
])

const ECL2ICRF = SMatrix{3, 3, Float64}([
    1.0 0.0                    0.0                  ;
    0.0 cosd(EARTH_OBLIQUITY) -sind(EARTH_OBLIQUITY);
    0.0 sind(EARTH_OBLIQUITY)  cosd(EARTH_OBLIQUITY);
])

