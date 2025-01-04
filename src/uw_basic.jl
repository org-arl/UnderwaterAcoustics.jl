import DSP: db2amp

export soundspeed, absorption, water_density, doppler, bubble_resonance
export reflection_coef, surface_reflection_coef, dBperλ, in_dBperλ

################################################################################
# common underwater acoustics utility functions

"""
    soundspeed(temperature=27, salinity=35, depth=0; γ=0, cₐ=340, ρᵣ=1000)

Compute sound speed in water in m/s, given:
- water `temperature` in °C
- `salinity` in ppt
- `depth` in meters
- void fraction (`γ`) in bubbly water
- sound speed in gas (`cₐ`) if `γ` > 0
- ratio of density of water to gas (`ρᵣ`) if `γ` > 0

Implementation based on Mackenzie (1981), Wood (1964) and Buckingham (1997).
"""
function soundspeed(temperature=27.0, salinity=35.0, depth=0.0; γ=0.0, cₐ=340.0, ρᵣ=1000.0)
  c = 1448.96 + 4.591*temperature - 5.304e-2*temperature^2 + 2.374e-4*temperature^3
  c += 1.340*(salinity-35) + 1.630e-2*depth + 1.675e-7*depth^2
  c += -1.025e-2*temperature*(salinity-35) - 7.139e-13*temperature*depth^3
  if γ > 0.0
    m = √ρᵣ
    c = 1.0/(1.0/c*√((γ*(c/cₐ)^2*m+(1.0-γ)/m)*(γ/m+(1-γ)*m)))
  end
  return c
end

"""
    absorption(frequency, distance=1000, salinity=35, temperature=27, depth=0, pH=8.1)

Compute volume acoustic absorption coefficient in water, given:
- `frequency` in Hz
- `distance` in meters
- `salinity` in ppm
- water `temperature` in °C
- `depth` in meters
- `pH` of water

The result is a unitless linear scale factor for sound pressure over the given `distance`. To get
absorption in terms of dB / m, set `distance = 1.0` and convert the result to decibels. For instance,
at a frequency of 100 kHz:

```julia-repl
julia> A = absorption(100e3, 1.0)
0.9959084838594522

julia> α = -20log10(A)
0.035611359656810865
```

Implementation based on the Francois and Garrison (1982) model.
"""
function absorption(frequency, distance=1000.0, salinity=35.0, temperature=27.0, depth=0.0, pH=8.1)
  f = frequency/1000.0
  d = distance/1000.0
  c = 1412.0 + 3.21*temperature + 1.19*salinity + 0.0167*depth
  A1 = 8.86/c * 10.0^(0.78*pH-5.0)
  P1 = 1.0
  f1 = 2.8*√(salinity/35.0) * 10.0^(4.0-1245.0/(temperature+273.0))
  A2 = 21.44*salinity/c*(1.0+0.025*temperature)
  P2 = 1.0 - 1.37e-4*depth + 6.2e-9*depth*depth
  f2 = 8.17 * 10^(8.0-1990.0/(temperature+273.0)) / (1.0+0.0018*(salinity-35.0))
  P3 = 1.0 - 3.83e-5*depth + 4.9e-10*depth*depth
  if temperature < 20.0
    A3 = 4.937e-4 - 2.59e-5*temperature + 9.11e-7*temperature^2 - 1.5e-8*temperature^3
  else
    A3 = 3.964e-4 - 1.146e-5*temperature + 1.45e-7*temperature^2 - 6.5e-10*temperature^3
  end
  a = A1*P1*f1*f*f/(f1*f1+f*f) + A2*P2*f2*f*f/(f2*f2+f*f) + A3*P3*f*f
  db2amp(-a*d)
end

"""
    water_density(temperature=27, salinity=35)

Compute density of water (kg/m^3), given `temperature` in °C and `salinity` in ppm.

Implementation based on Fofonoff (1985 - IES 80).
"""
function water_density(temperature=27.0, salinity=35.0)
  t = temperature
  A = 1.001685e-04 + t * (-1.120083e-06 + t * 6.536332e-09)
  A = 999.842594 + t * (6.793952e-02 + t * (-9.095290e-03 + t * A))
  B = 7.6438e-05 + t * (-8.2467e-07 + t * 5.3875e-09)
  B = 0.824493 + t * (-4.0899e-03 + t * B)
  C = -5.72466e-03 + t * (1.0227e-04 - t * 1.6546e-06)
  D = 4.8314e-04
  A + salinity * (B + C*√(salinity) + D*salinity)
end

"""
    reflection_coef(θ, ρᵣ, cᵣ, δ=0.0)

Compute complex reflection coefficient at a fluid-fluid boundary, given:
- angle of incidence `θ` (angle to the surface normal)
- relative density of the reflecting medium to incidence medium `ρᵣ`
- relative sound speed of the reflecting medium to incidence medium `cᵣ`
- dimensionless absorption coefficient `δ`

Implementation based on Brekhovskikh & Lysanov. Dimensionless absorption
coefficient based on APL-UW Technical Report 9407.
"""
function reflection_coef(θ, ρᵣ, cᵣ, δ=0.0)
  n = Complex(1.0, δ) / cᵣ
  t1 = ρᵣ * cos(θ)
  t2 = n*n - sin(θ)^2
  t3 = √abs(t2) * cis(angle(t2)/2)   # ForwardDiff friendly complex √
  (t1 - t3) / (t1 + t3)
end

"""
    dBperλ(x)

Compute dimensionless absorption coefficient `δ` from dB/λ.
Implementation based on APL-UW TR 9407 (1994), IV-9 equation (4).
"""
dBperλ(x) = x / (40π / log(10))

"""
    in_dBperλ(δ)

Compute dB/λ from dimensionless absorption coefficient `δ`.
Implementation based on APL-UW TR 9407 (1994), IV-9 equation (4).
"""
in_dBperλ(δ) = δ * 40π / log(10)

"""
    surface_reflection_coef(windspeed, frequency, θ)


Compute surface reflection coefficient, given:
- `windspeed` in m/s
- `frequency` in Hz
- angle of incidence `θ` (angle to the surface normal)

Implementation based on the APL-UW Technical Report 9407 II-21.
"""
function surface_reflection_coef(windspeed, frequency, θ)
  β = π/2 - θ
  f = frequency/1000.0
  if windspeed >= 6.0
    a = 1.26e-3/sin(β) * windspeed^1.57 * f^0.85
  else
    a = 1.26e-3/sin(β) * 6^1.57 * f^0.85 * exp(1.2*(windspeed-6.0))
  end
  db2amp(-a)
end

"""
    doppler(speed, frequency)
    doppler(speed, frequency, soundspeed)

Compute Doppler frequency, given relative speed between transmitter and
receiver in m/s. `soundspeed` is the nominal sound speed in water.
"""
doppler(speed, frequency, soundspeed=soundspeed()) = (1 + speed / soundspeed) * frequency

"""
    bubble_resonance(radius, depth=0, γ=1.4, p₀=1.013e5, ρ=1022.476)

Compute resonance frequency of a freely oscillating has bubble in water, given:
- bubble `radius` in meters
- `depth` of bubble in water in meters
- gas ratio of specific heats 'γ', default: 1.4 (for air)
- atmospheric pressure 'p₀', default: 1.013e5
- density of water 'ρ' in kg/m³, default: 1022.476

This ignores surface-tension, thermal, viscous and acoustic damping effects, and the pressure-volume relationship is taken to be adiabatic.
Implementation based on Medwin & Clay (1998).
"""
function bubble_resonance(radius, depth=0.0, γ=1.4, p₀=1.013e5, ρ=1022.476)
  g = 9.80665 # acceleration due to gravity
  pₐ = p₀ + ρ*g*depth
  1 / (2π * radius) * √(3γ * pₐ/ρ)
end
