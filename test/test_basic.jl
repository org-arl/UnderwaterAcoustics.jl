using TestItems

@testsnippet BasicSetup begin
  in_dB(x) = 20 * log10(x)
  dB(x) = 10^(x/20)
  const ° = u"°"
end

@testitem "soundspeed" begin
  @test @inferred(soundspeed()) ≈ 1539.0 atol=0.1
  @test @inferred(soundspeed(; γ=1e-5)) ≈ 1402.1 atol=0.1
  @test @inferred(soundspeed(; γ=1.0)) ≈ 340.0
  @test @inferred(soundspeed(27u"°C", 35u"ppt", 0u"m")) ≈ 1539.0 atol=0.1
  @test @inferred(soundspeed(; γ=1e-5, cₐ=340u"m/s")) ≈ 1402.1 atol=0.1
end

@testitem "shearspeed" begin
  @test @inferred(shearspeed(1520)) ≈ 146.68
  @test @inferred(shearspeed(1560)) ≈ 288.72
  @test @inferred(shearspeed(2000)) ≈ 599.0
  @test @inferred(shearspeed(2200)) ≈ 754.0
end

@testitem "absorption" setup=[BasicSetup] begin
  @test @inferred(in_dB(absorption(10000, 1000.0, 35.0, 15.0))) ≈ -1.0 atol=0.5
  @test @inferred(in_dB(absorption(50000))) ≈ -11.0 atol=0.5
  @test @inferred(in_dB(absorption(100000))) ≈ -36.0 atol=0.5
  @test @inferred(in_dB(absorption(100000, 1000.0, 38.5, 14.0, 0.0))) ≈ -40.0 atol=0.5
  @test @inferred(in_dB(absorption(100000, 1000.0, 38.5, 14.0, 2000.0))) ≈ -30.0 atol=0.5
  @test @inferred(in_dB(absorption(100000, 1000.0, 38.5, 14.0, 6000.0))) ≈ -16.0 atol=0.5
  @test @inferred(in_dB(absorption(100u"kHz", 1u"km", 38.5u"ppt", 14u"°C", 6u"km"))) ≈ -16.0 atol=0.5
end

@testitem "water_density" setup=[BasicSetup] begin
  @test @inferred(water_density()) ≈ 1022.7 atol=0.1
  @test @inferred(water_density(27u"°C", 35u"ppt")) ≈ 1022.7 atol=0.1
end

@testitem "doppler" begin
  @test @inferred(doppler(0.0, 50000.0)) == 50000.0
  @test @inferred(doppler(10.0, 50000.0)) ≈ 50325 atol=0.5
  @test @inferred(doppler(-10.0, 50000.0)) ≈ 49675 atol=0.5
  @test @inferred(doppler(-10u"m/s", 50u"kHz")) ≈ 49675 atol=0.5
end

@testitem "bubble_resonance" begin
  @test @inferred(bubble_resonance(100e-6)) ≈ 32461.7 atol=0.1
  @test @inferred(bubble_resonance(32e-6)) ≈ 101442.8 atol=0.1
  @test @inferred(bubble_resonance(100e-6, 10.0)) ≈ 45793.7 atol=0.1
  @test @inferred(bubble_resonance(100u"µm", 10u"m"; p₀=101.3u"kPa", ρ=1022.72u"kg/m^3", g=9.80665u"m/s^2")) ≈ 45793.7 atol=0.1
end

@testitem "reflection_coef" setup=[BasicSetup] begin
  @test @inferred(reflection_coef(0.0, 1200.0/1023.0, 1600.0/1540.0, 0.0)) isa Complex
  @test @inferred(reflection_coef(0.0, 1200.0/1023.0, 1600.0/1540.0, 0.0)) ≈ 0.0986 atol=0.0001
  @test @inferred(in_dB(abs(reflection_coef(0.0, 2.5, 2.5, 0.01374)))) ≈ -2.8 atol=0.1
  @test @inferred(in_dB(abs(reflection_coef(0.0, 2.492, 1.3370, 0.01705)))) ≈ -5 atol=0.5
  @test @inferred(in_dB(abs(reflection_coef(0.0, 1.195, 1.0179, 0.02158)))) ≈ -20 atol=0.5
  @test @inferred(in_dB(abs(reflection_coef(1.22, 1.149, 0.9873, 0.00386)))) ≈ -32 atol=0.5
  @test @inferred(in_dB(abs(reflection_coef(69°, 1.149, 0.9873, 0.00386)))) ≈ -32 atol=0.5
end

@testitem "surface_reflection_coef" setup=[BasicSetup] begin
  @test @inferred(surface_reflection_coef(15.0, 20000.0, 80°)) isa Complex
  @test @inferred(surface_reflection_coef(15.0, 20000.0, 80°)) ≈ -dB(-6.5) atol=0.1
  @test @inferred(surface_reflection_coef(10.0, 20000.0, 80°)) ≈ -dB(-3.4) atol=0.1
  @test @inferred(surface_reflection_coef(5.0, 20000.0, 80°)) ≈ -dB(-0.5) atol=0.1
  @test @inferred(surface_reflection_coef(5u"m/s", 20u"kHz", 1.3962634)) ≈ -dB(-0.5) atol=0.1
end

@testitem "dBperλ" begin
  @test @inferred(dBperλ(0.0)) == 0.0
  @test @inferred(dBperλ(40π/log(10))) ≈ 1.0
  @test @inferred(in_dBperλ(0.0)) == 0.0
  @test @inferred(in_dBperλ(1.0)) ≈ 40π/log(10)
  @test @inferred(in_dBperλ(dBperλ(10.0))) ≈ 10.0
end
