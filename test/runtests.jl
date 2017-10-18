#module PowerFlowTests
include("../src/PowerFlowAnalysis.jl")

using Base.Test
using JuMP
using Ipopt
using PowerModels
using PowerFlowAnalysis

case6ww = PowerModels.parse_file("data/case6ww.m")
pm = ACPPowerModel(case6ww)
sol = run_ac_pf(case6ww, IpoptSolver())

net = Network(pm)
sys = PowerFlowSystem(pm, sol)
lf = LurieForm(sys)

known_pairs = [(1,4), (2,5), (3,6)] # lines in the orignal file
@testset "Network topology" for pair in known_pairs
  f, t = pair
  @test (f,t) in net.pairs
  pidx = net.pair_idx[(f,t)]
  @test net.from[pidx] == net.bus_idx[f]
  @test net.to[pidx]   == net.bus_idx[t]
end

@testset "Admittance matrix" begin
  Ymp = Dict(
      (5,5) => 5.293290152818541 - 14.137764902337823im,
      (3,6) => -1.923076923076923 + 9.615384615384615im,
      (1,3) => 0.0im
    ) # Several values extracted from MATPOWER
  for ((f,t), y) in Ymp
    @test y ≈ sys.Y[net.bus_idx[f], net.bus_idx[t]]
  end
end

p, q = calc_pq(sys.Y, sys.v⁰, sys.ϑ⁰)
@testset "Power flow solutions" begin
  @test p[sys.pq] ≈ sys.p⁰[sys.pq]
  @test p[sys.pv] ≈ sys.p⁰[sys.pv]
  @test q[sys.pq] ≈ sys.q⁰[sys.pq]
end

@testset "Lurie form" begin
  v¹ = (1 + 0.1*randn(sys.nb)).*sys.v⁰
  ϑ¹ = (1 + 0.1*randn(sys.nb)).*sys.ϑ⁰
  f¹ = calc_f(v¹, ϑ¹, lf)
  s¹ = lf.S*f¹
  p¹, q¹ = calc_pq(sys.Y, v¹, ϑ¹)

  @test p¹[sys.pv] ≈ s¹[1:sys.npv]
  @test p¹[sys.pq] ≈ s¹[(sys.npv + 1):(sys.npv + sys.npq)]
  @test q¹[sys.pq] ≈ s¹[(sys.npv + sys.npq + 1):(sys.npv + 2*sys.npq)]
end

@testset "Jacobian ∂f/∂x" begin
  dfdx = PowerFlowAnalysis.calc_dfdx(lf)
  dx = 1e-4*randn(size(lf.x⁰))
  v¹,ϑ¹ = PowerFlowAnalysis.calc_vϑ(lf.x⁰ + dx, lf)
  df = PowerFlowAnalysis.calc_f(v¹, ϑ¹, lf) - lf.f⁰
  @test isapprox(df, dfdx*dx; rtol = 1e-3)
end
