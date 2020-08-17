export soundspeed, absorption, density, reflectioncoef, surfaceloss, doppler, bubbleresonance

function soundspeed(temperature=27.0, salinity=35.0, depth=10.0; ν=0.0, cgas=340.0, reldensity=1000.0)
  # based on Mackenzie (1981)
  c = 1448.96 + 4.591*temperature - 5.304e-2*temperature^2 + 2.374e-4*temperature^3
  c += 1.340*(salinity-35) + 1.630e-2*depth + 1.675e-7*depth^2
  c += -1.025e-2*temperature*(salinity-35) - 7.139e-13*temperature*depth^3
  if ν > 0.0
    # based on Wood (1964) or Buckingham (1997)
    m = √reldensity
    c = 1.0/(1.0/c*√((ν*(c/cgas)^2*m+(1.0-ν)/m)*(ν/m+(1-ν)*m)))
  end
  return c
end

function absorption(frequency, distance=1000.0, salinity=35.0, temperature=27.0, depth=10.0, pH=8.1)
  # based on Francois-Garrison model
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
  10^(-a*d/20.0)
end

function waterdensity(temperature=27, salinity=35)
  # based on Fofonoff (1985 - IES 80)
  t = temperature
  A = 1.001685e-04 + t * (-1.120083e-06 + t * 6.536332e-09)
  A = 999.842594 + t * (6.793952e-02 + t * (-9.095290e-03 + t * A))
  B = 7.6438e-05 + t * (-8.2467e-07 + t * 5.3875e-09)
  B = 0.824493 + t * (-4.0899e-03 + t * B)
  C = -5.72466e-03 + t * (1.0227e-04 - t * 1.6546e-06)
  D = 4.8314e-04
  A + salinity * (B + C*√(salinity) + D*salinity)
end

function reflectioncoef(θ, ρᵣ, cᵣ, δ=0.0)
  # based on Brekhovskikh & Lysanov
  n = Complex(1.0, δ) / cᵣ
  t1 = ρᵣ * cos(θ)
  t2 = √(n^2 - sin(θ)^2)
  (t1-t2)/(t1+t2)
end

function surfaceloss(windspeed, frequency, θ)
  # based on APL model (1994)
  β = π/2 - θ
  if windspeed >= 6.0
    a = 1.26e-3/sin(β) * windspeed^1.57 * frequency^0.85
  else
    a = 1.26e-3/sin(β) * 6^1.57 * frequency^0.85 * exp(1.2*(windspeed-6.0))
  end
  db2amp(-a)
end

doppler(speed, frequency, c=soundspeed()) = (1.0+speed/c)*frequency
bubbleresonance(radius, depth=0) = 3.25/radius * √(1+0.1*depth)     # Medwin & Clay (1998)
