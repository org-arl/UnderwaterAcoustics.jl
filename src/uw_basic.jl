export soundspeed, absorption, waterdensity, reflectioncoef, surfaceloss, doppler, bubbleresonance

"""
$(SIGNATURES)
Compute sound speed in water in m/s, given:
- water `temperature` in °C
- `salinity` in ppm
- `depth` in meters
- void fraction (`voidfrac`) in bubbly water
- sound speed in gas (`cgas`), if `voidfrac` > 0
- ratio of density of water to gas (`reldensity`), if `voidfrac` > 0

Implementation based on Mackenzie (1981), Wood (1964) and Buckingham (1997).
"""
function soundspeed(temperature=27.0, salinity=35.0, depth=10.0; voidfrac=0.0, cgas=340.0, reldensity=1000.0)
  c = 1448.96 + 4.591*temperature - 5.304e-2*temperature^2 + 2.374e-4*temperature^3
  c += 1.340*(salinity-35) + 1.630e-2*depth + 1.675e-7*depth^2
  c += -1.025e-2*temperature*(salinity-35) - 7.139e-13*temperature*depth^3
  if voidfrac > 0.0
    m = √reldensity
    c = 1.0/(1.0/c*√((voidfrac*(c/cgas)^2*m+(1.0-voidfrac)/m)*(voidfrac/m+(1-voidfrac)*m)))
  end
  return c
end

"""
$(SIGNATURES)
Compute volume acoustic absorption coefficient in water, given:
- `frequency` in Hz
- `distance` in meters
- `salinity` in ppm
- water `temperature` in °C
- `depth` in meters
- `pH` of water

Implementation based on the Francois-Garrison model.
"""
function absorption(frequency, distance=1000.0, salinity=35.0, temperature=27.0, depth=10.0, pH=8.1)
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
$(SIGNATURES)
Compute density of water (kg/m^3), given `temperature` in °C and `salinity` in ppm.

Implementation based on Fofonoff (1985 - IES 80).
"""
function waterdensity(temperature=27, salinity=35)
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
$(SIGNATURES)
Compute complex reflection coefficient at a fluid-fluid boundary, given:
- angle of incidence `θ` (angle to the surface normal)
- relative density of the reflecting medium to incidence medium `ρᵣ`
- relative sound speed of the reflecting medium to incidence medium `cᵣ`
- dimensionless absorption coefficient `δ`

Implementation based on Brekhovskikh & Lysanov. Dimensionless absorption
coefficient based on APL-UW Technical Report 9407.
"""
function reflectioncoef(θ, ρᵣ, cᵣ, δ=0.0)
  n = Complex(1.0, δ) / cᵣ
  t1 = ρᵣ * cos(θ)
  t2 = n*n - sin(θ)^2
  t3 = √abs(t2) * cis(angle(t2)/2)   # ForwardDiff friendly complex √
  (t1 - t3) / (t1 + t3)
end

"""
$(SIGNATURES)
Compute surface reflection coefficient, given:
- `windspeed` in m/s
- `frequency` in Hz
- angle of incidence `θ` (angle to the surface normal)

Implementation based on the APL-UW Technical Report 9407 II-21.
"""
function surfaceloss(windspeed, frequency, θ)
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
$(SIGNATURES)
Compute Doppler frequency, given relative speed between transmitter and
receiver in m/s. `c` is the nominal sound speed in water.
"""
doppler(speed, frequency, c=soundspeed()) = (1.0+speed/c)*frequency

"""
$(SIGNATURES)
Compute resonance frequency of a freely oscillating has bubble in water, given:
- bubble `radius` in meters
- `depth` of bubble in water in meters
- gas ratio of specific heats 'γ', default: 1.4 (for air)
- atmospheric pressure 'p0' in Pa, default: 1.013e5
- density of water 'ρ' in kg/m³, default: 1022.476

This ignores surface-tension, thermal, viscous and acoustic damping effects, and the pressure-volume relationship is taken to be adiabatic.
Implementation based on Medwin & Clay (1998).
"""
function bubbleresonance(radius::Real, depth::Real=0.0, γ::Real=1.4, p0::Real=1.013e5, ρ::Real=1022.476)
  g = 9.80665 #acceleration due to gravity
  pₐ = p0 + ρ*g*depth
  1 / (2π * radius) * √(3γ * pₐ/ρ)
end

"""
$(SIGNATURES)
Compute bubble resonant frequency and damping constants. Uses Prosperetti (1977) equations to compute damping coefficients and natural frequency.
(Thermal effects and damping mechanisms in the forced radial oscillations of gas bubbles in liquids).
  
Parameters:
- bubble 'radius' in meters
- frequency 'f' in Hz
- water depth 'z' in meters. default: 0
- bubble gas and water temperature 'T' in Kelvin. default: 293.15
- atmospheric pressure 'p0' in Pa. default: 101.3e5
- 'gas' number. If 2: Oxygen, 3: Nitrogen, 4: Methane, any other number: Air. default: 1
- Soundspeed 'c' in m/s. default: 1500
- seawater density 'ρ' in kg m^-3. default: 1022.476

Returns a 1x4 Array{Float64}
- ω0   : bubble natural frequency in radians/s
- b_th : thermal damping coefficient in /s
- b_ac : radiation damping coefficient in /s
- b_vs : viscous damping coefficient in /s
"""  
using LinearAlgebra
function bubble_resonance_and_damping(radius::Real, f::Real; z::Real=0., T::Real=293.15, p0=1.013e5, gas::Integer = 1, c::Real=1500., ρ::Real=1022.476)
    ## Set up the physical constants
    g::Float64 = 9.80665;        # acceleration due to gravity
    Rg::Float64 = 8.314472;      # universal gas constant in Joules K^-1 mol^-1
    # Water Properties. These are assumed independent of pressure
    # Seawater density.
    # http://www.es.flinders.edu.au/~mattom/Utilities/density.html
    # T = 20 Celsius, salinity = 32 ppt, pressure = 0
    Kw::Float64 = 0.58;          # thermal conductivity water W/m/K
    Dw::Float64 = 1.45e-7;       # thermal diffusivity of water in m^2 s^-1
    μ::Float64 = 1e-3;          # viscosity of water in Pa.s
    σ::Float64 = 0.072;        # water surface tension in N m
  
    # Ambient pressure and internal, equilibrium bubble pressure at depth
    pamb::Float64 = p0 + z*ρ*g; # ambient water pressure - pamb is denoted as p0 in prosperetti 76a
    pineq::Float64 = pamb .+ 2σ./radius; #pg 3. below eq 6
    ω::Float64 = 2π*f;
    k::Float64 = ω/c;
    # Set up the gas values
    # The thermal conductivity and ratio of specific heats came from:
    # http://www.engineeringtoolbox.com/oxygen-d_978.html
    # The thermal conductivity came from:
    # http://www.webelements.com/webelements/elements/text/O/heat.html
    # M::Float64
    # Kg::Float64
    # Cp::Float64
    # γ::Float64
    if gas == 2 # Oxygen
      M = 0.0320;     # Molecular weight Oxygen Kg / mol
      Kg = 0.02658;    # Thermal conductivity Oxygen W /m /K
      Cp = 915 + (918-915)/(300-275)*(290-275);     # J/Kg/K (290 K)
      γ = 1.4;     # ratio of specific heats
    elseif gas == 3 # Nitrogen
      M =  0.02802;  # Molecular weight Nitrogen Kg / mol
      Kg = 0.02583;  # Thermal conductivity Nitrogen W /m /K
      Cp = 1039 + (1040-1039)/(300-275)*(290-275); # J/Kg/K (290 K)
      γ = 1.4;   # ratio of specific hetas
    elseif gas == 4  #Methane
      M = 0.016;  # Molecular weight methane Kg / mol
      Kg = 0.0332; # Thermal conductivity methane W /m /K (290 K)
      Cp = 2191 + (2226-2191)/(300-275)*(290-275);  # J/Kg/K (290 K)
      γ = 1.31;   # ratio of specific hetas
    else  #Air
      # Specific Heat Ratio for air.  http://www.efunda.com/Materials/common_matl/show_gas.cfm?MatlName=AiradiusC
      # Thermal conductivity of air. This is assumed to be independent of
      # pressure. http://home.worldonline.dk/jsrsw/Tcondvspressure.html
      M = 0.02896;       # air molecular weight in kg mol^-1 
      Kg = 0.0254;        # thermal conductivity air W/m/K
      Cp = 1.00e3;        # specific heat capacity of air in J Kg^-1 K^-1
      γ = 1.40;       # ratio of specific heats
    end
    
    # Compute bubble gas density from the internal, equilibrium pressure
    ρg::Float64 = M*pineq/(Rg*T);
    # Compute thermal diffusivity of gas. Note Prosperetti (1976a, p.19) has a 
    # typo: Cv,g should be Cp,g as used beloω.
    Dg::Float64 = Kg/(ρg*Cp);
    ## Losses versus frequency for fixed radius and depth.
    G1::Float64 = M*Dg*ω/(γ*Rg*T); #pg 19
    G2::Float64 = ω.*(radius.^2)./Dg; #pg 19. essentialy the square of the ratio between the bubble tabus and the thermal penetration depth
    G3::Float64 = ω.*(radius.^2)./Dw; #pg 19
    F::Number = 1 .+ (1+1im)*.√(0.5*G3);
    kratio::Float64 = Kw/Kg;
  
    β1::Number = .√(0.5γ*G2.* (1im .- G1 .+ .√((1im .- G1).^2 .+ 4im*G1/γ)));
    β2::Number = .√(0.5γ*G2.* (1im .- G1 .- .√((1im .- G1).^2 .+ 4im*G1/γ)));
    λ1::Number = β1.*coth.(β1) .- 1;
    λ2::Number = β2.*coth.(β2) .- 1;
    Γ1::Number = 1im .+ G1 .+ .√((1im .- G1).^2 .+ 4im*G1/γ);
    Γ2::Number = 1im .+ G1 .- .√((1im .- G1).^2 .+ 4im*G1/γ);
    ψ::Number = (kratio*F.*(Γ2 .- Γ1) .+ λ2.*Γ2 .- λ1.*Γ1)./(kratio.*F.*(λ2.*Γ1 .- λ1.*Γ2)- λ1.*λ2.*(Γ2.-Γ1)); 
  
    μth::Float64 = 1/4*(ω.*ρg.*radius.^2).*imag.(ψ);#effective thermal viscosity. 
    κ::Float64 = 1/3*((ω.^2).*ρg.*radius.^2 ./pineq).*real.(ψ); #polytropic exponent. 
    b_th::Float64 = 2μth./(ρ.*radius.^2); #thermal damping
    b_ac::Float64 = 1/2 .* ω.* (ω.*radius/c)./(1 .+ (ω.*radius/c).^2); #acoustic damping
    # b_vs = repeat([2μ./(ρ*radius.^2)], outer = [1,length(ω)]); #viscous damping
    b_vs::Float64 = 2μ./(ρ.*radius.^2); #viscous damping
    ω0::Float64 = .√(3κ .* pineq./(ρ.*radius.^2) .- 2σ./(ρ.*(radius.^3)))# + (k.^2).*(radius.^2)./(1 .+ (k.^2).*(radius.^2)).*(ω.^2) );
    return [ω0,b_th, b_ac, b_vs]
end

""" 
$(SIGNATURES)
Compute Scattering coefficient due to bubble suspended in water. 
Damping effects are included, and no assumptions about spherical nature of bubbles made.

Parameters:
- bubble 'radius' in meters
- frequency vector 'f' in Hz 
- bubble 'depth' default: 0 m
- water temperature 'T' in Kelvin. default 293.15K
- atmospheric pressure 'p0' in Pa. default: 1.013e5 Pa
- gas number. If 2: Oxygen, 3: Nitrogen, 4: Methane, any other number: Air. Default is 1.
- soundspeed 'c' in m/s. default: 1500
- seawater density 'ρ' in kg m^-3. default: 1022.476
"""
function scattering_coeff(radius::Real, f::Real; depth::Real=0, T::Real=293.15, p0::Real=1.013e5, gas::Integer=1, c::Real=1500, ρ::Real=1022.476)
  ω = 2π*f

  bubbleparams = bubble_resonance_and_damping(radius, f, z=depth, T=T, p0=p0, gas=gas, c=c)
  ω0 = bubbleparams[1]
  b_th = bubbleparams[2]
  b_ac = bubbleparams[3]
  b_vs = bubbleparams[4]
  β = b_th .+ b_ac .+ b_vs #sum of all 3 types of damping - thermal, viscous and acoustic
  denom = (ω0./ω).^2 .- 1 .+ 2im*β.*(1 ./ω)
  gscat = radius./denom
  return gscat
end
