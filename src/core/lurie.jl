type LurieFormAd
  sys::PowerFlowSystem

  x⁰::Vector
  f⁰::Vector
  s⁰::Vector

  S::AbstractMatrix
  J⁰::AbstractMatrix
  J⁻¹::AbstractMatrix

  SQ::AbstractMatrix
  SQJ::AbstractMatrix

  function LurieFormAd(sys::PowerFlowSystem)
    @import_fields sys

    x⁰ = calc_x(v⁰, ϑ⁰, sys)
    f⁰ = calc_f(v⁰, ϑ⁰, sys)

    dfdx⁰ = calc_dfdx(v⁰, ϑ⁰, sys)
    S = calc_S(sys)

    J⁰ = S*dfdx⁰
    s⁰ = S*f⁰
    J⁻¹ = inv(full(J⁰))

    SQ = calc_S_Q(sys)
    SQJ = SQ*dfdx⁰

    new(sys, x⁰, f⁰, s⁰, S, J⁰, J⁻¹, SQ, SQJ)
  end
end

LurieFormAd(pm::GenericPowerModel) = LurieFormAd(PowerFlowSystem(pm))

calc_x(v::Vector, ϑ::Vector, lf::LurieFormAd) = calc_x(v, ϑ, lf.sys)
function calc_x(v::Vector, ϑ::Vector, sys::PowerFlowSystem)
  @import_fields sys
  return [ϑ[pv]; ϑ[pq]; log(v[pq])]
end

calc_ρᵦϑ(x::Vector, lf::LurieFormAd) = calc_ρᵦϑ(x, lf.sys) #rho of nodes
function calc_ρᵦϑ(x::Vector, sys::PowerFlowSystem)
  @import_fields sys
  v = copy(v⁰); ϑ = copy(ϑ⁰); ρᵦ = copy(v⁰);
  ϑ[pv] = x[1:npv]
  ϑ[pq] = x[(npv + 1):(npv + npq)]
  v[pq] = x[(npv + npq + 1):(npv + 2*npq)]
  ρᵦ[pq] = log(v[pq])
  return (ρᵦ, ϑ)
end

calc_f(v::Vector, ϑ::Vector, lf::LurieFormAd) = calc_f(v, ϑ, lf.sys)
function calc_f(v::Vector, ϑ::Vector, sys::PowerFlowSystem)
  @import_fields sys
  ρₑ = log(v[from]./v[to]) - log(v⁰[from]./v⁰[to]) #rho of edges
  Θ = (ϑ[from] - ϑ[to]) - (ϑ⁰[from] - ϑ⁰[to])
  [cosh(ρₑ).*cos(Θ) ; sinh(ρₑ).*cos(Θ); cosh(ρₑ).*sin(Θ); sinh(ρₑ).*sin(Θ)]
end

range_coshcos(lf::LurieFormAd) = range_coshcos(lf.sys.np)
range_coshcos(sys::PowerFlowSystem) = range_coshcos(sys.np)
range_coshcos(np::Integer) = collect(1:np)

range_sinhcos(lf::LurieFormAd) = range_sinhcos(lf.sys.np)
range_sinhcos(sys::PowerFlowSystem) = range_sinhcos(sys.np)
range_sinhcos(np::Integer) = collect((np+1):(2*np))

range_coshsin(lf::LurieFormAd) = range_coshsin(lf.sys.np)
range_coshsin(sys::PowerFlowSystem) = range_coshsin(sys.np)
range_coshsin(np::Integer) = collect((2*np+1):(3*np))

range_sinhsin(lf::LurieFormAd) = range_sinhsin(lf.sys.np)
range_sinhsin(sys::PowerFlowSystem) = range_sinhsin(sys.np)
range_sinhsin(np::Integer) = collect((3*np+1):(4*np))

calc_dfdx(lf::LurieFormAd) = calc_dfdx(lf.sys)
calc_dfdx(sys::PowerFlowSystem) = calc_dfdx(sys.v⁰, sys.ϑ⁰, sys)

function calc_dfdx(v::Vector, ϑ::Vector, sys::PowerFlowSystem)
  @import_fields sys

  f = calc_f(v, ϑ, sys)
  coshcos = f[range_coshcos(sys)]
  sinhcos = f[range_sinhcos(sys)]
  coshsin = f[range_coshsin(sys)]
  sinhsin = f[range_sinhsin(sys)]
  frto = [from;to]
  pvpq = [pv;pq]
  np2  = [1:np;1:np]

  ∂chcs╱∂ρ = sparse(np2, frto, [sinhcos; -sinhcos], np, nb)
  ∂chcs╱∂ϑ = sparse(np2, frto, [-coshsin; coshsin], np, nb)

  ∂shcs╱∂ρ = sparse(np2, frto, [coshcos; -coshcos], np, nb)
  ∂shcs╱∂ϑ = sparse(np2, frto, [-sinhsin; sinhsin], np, nb)

  ∂chsn╱∂ρ = sparse(np2, frto, [sinhsin; -sinhsin], np, nb)
  ∂chsn╱∂ϑ = sparse(np2, frto, [coshcos; -coshcos], np, nb)

  ∂shsn╱∂ρ = sparse(np2, frto, [coshsin; -coshsin], np, nb)
  ∂shsn╱∂ϑ = sparse(np2, frto, [sinhcos; -sinhcos], np, nb)

  return [
    ∂chcs╱∂ϑ[:,pvpq]    ∂chcs╱∂ρ[:,pq]
    ∂shcs╱∂ϑ[:,pvpq]    ∂shcs╱∂ρ[:,pq]
    ∂chsn╱∂ϑ[:,pvpq]    ∂chsn╱∂ρ[:,pq]
    ∂shsn╱∂ϑ[:,pvpq]    ∂shsn╱∂ρ[:,pq]
         ]
end

# There is some bug with the Base.real implementation
spreal(A::SparseMatrixCSC) = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, real(A.nzval))
spimag(A::SparseMatrixCSC) = SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval, imag(A.nzval))

#
function calc_S(sys::PowerFlowSystem)
  @import_fields sys
  ρ⁰ = log(v⁰[from] ./ v⁰[to])
  Θ⁰ = ϑ⁰[from] - ϑ⁰[to]
  Ŷᵗ = Yᵗ * spdiagm(exp(ρ⁰ + Θ⁰*1im), (0))
  Ŷᶠ = Yᶠ * spdiagm(exp(-ρ⁰ - Θ⁰*1im), (0))
  Y⁺ = Ŷᶠ + Ŷᵗ
  Y⁻ = Ŷᶠ - Ŷᵗ

  S = [
        spreal(Y⁺[pv,:])  spreal(Y⁻[pv,:])  -spimag(Y⁻[pv,:])  -spimag(Y⁺[pv,:])
        spreal(Y⁺[pq,:])  spreal(Y⁻[pq,:])  -spimag(Y⁻[pq,:])  -spimag(Y⁺[pq,:])
        -spimag(Y⁺[pq,:])  -spimag(Y⁻[pq,:])  -spreal(Y⁻[pq,:])  -spimag(Y⁺[pq,:])
      ]
  return S
end

function calc_S_Q(sys::PowerFlowSystem)
  @import_fields sys
  ρ⁰ = log(v⁰[from] ./ v⁰[to])
  Θ⁰ = ϑ⁰[from] - ϑ⁰[to]
  Ŷᵗ = Yᵗ * spdiagm(exp(ρ⁰ + Θ⁰*1im), (0))
  Ŷᶠ = Yᶠ * spdiagm(exp(-ρ⁰ - Θ⁰*1im), (0))
  Y⁺ = Ŷᶠ + Ŷᵗ
  Y⁻ = Ŷᶠ - Ŷᵗ

  SQ = [spimag(Y⁺[pv,:])  spimag(Y⁻[pv,:])  -spreal(Y⁻[pv,:])  spimag(Y⁺[pv,:])]
  return SQ
end

export LurieFormAd, calc_x, calc_vϑ, calc_f, calc_dfdx
export range_coshcos, range_coshsin, range_sinhcos, range_sinhsin
