typealias IndexList Vector{Int}
typealias ComplexMatrix AbstractMatrix{Complex}
#typealias Vector Vector

type PowerFlowSystem
  buses::Vector # Vector of bus names/indices

  pq::IndexList # pq bus indices
  pv::IndexList # pv bus indices
  sl::IndexList # slack bus indices

  from::IndexList # list of pair "from" bus indices
  to::IndexList   # list of pair "to"  bus indices

  Yᵈ::ComplexMatrix # diagonal of Y matrix
  Yᶠ::ComplexMatrix # admittance of edges exiting the bus
  Yᵗ::ComplexMatrix # admittance of edges entering the bus
  Y::ComplexMatrix  # admittance matrix

  v⁰::Vector # base solution voltage magnitudes
  ϑ⁰::Vector # base solution voltage angles

  p⁰::Vector # base solution active powers
  q⁰::Vector # base solution reactive powers

  nb::Integer # total number of buses
  np::Integer # total number of pairs

  npq::Integer # total number of PQ buses
  npv::Integer # total number of PV buses
  nsl::Integer # total number of Slack buses

  function PowerFlowSystem(buses, pq, pv, slack, from, to, Yᵈ, Yᶠ, Yᵗ, Y, v⁰, ϑ⁰ , p⁰, q⁰)
    nb = length(buses)
    np = length(from)
    npq = length(pq)
    npv = length(pv)
    nsl = length(slack)
    new(buses, pq, pv, slack, from, to, Yᵈ, Yᶠ, Yᵗ, Y, v⁰, ϑ⁰, p⁰, q⁰, nb, np, npq, npv, nsl)
  end
end

function PowerFlowSystem(pm::GenericPowerModel)
  net = Network(pm)
  bt = get_bus_types(pm, net)
  YYY = calc_YYY(pm, net)
  Y = calc_Y(YYY..., net)
  vϑ = get_vϑ(pm, net)
  pq = get_pq(pm, net)
  return PowerFlowSystem(net.buses, bt..., net.from, net.to, YYY..., Y, vϑ..., pq...)
end

function PowerFlowSystem(pm::GenericPowerModel, sol::Dict)
  net = Network(pm)
  bt = get_bus_types(pm, net)
  YYY = calc_YYY(pm, net)
  Y = calc_Y(YYY..., net)
  vϑ = get_sol_vϑ(sol, net)
  pq = get_pq(pm, net)
  return PowerFlowSystem(net.buses, bt..., net.from, net.to, YYY..., Y, vϑ..., pq...)
end

function calc_pq(Y::AbstractMatrix, v::Vector, ϑ::Vector)
  V = v.*exp(im*ϑ)
  return reim(V.*conj(Y*V))
end

function calc_Y(Yᵈ::AbstractMatrix, Yᶠ::AbstractMatrix, Yᵗ::AbstractMatrix, net::Network)
  @import_fields net
  Cf = sparse(1:np,   to, 1.0 + 0im, np, nb)
  Ct = sparse(1:np, from, 1.0 + 0im, np, nb)
  return Yᵈ + Yᶠ*Cf + Yᵗ*Ct
end

export PowerFlowSystem, calc_pq, calc_Y
