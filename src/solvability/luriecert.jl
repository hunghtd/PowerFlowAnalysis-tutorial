using MATLAB

type BaseLurieAdCertificate
  lf::LurieFormAd

  B::AbstractMatrix
  C::AbstractMatrix
  R::AbstractMatrix

  D::AbstractMatrix
  E::AbstractMatrix

  function BaseLurieAdCertificate(lf::LurieFormAd, R)
    A = create_A(lf)
    B, C = create_BC(lf, A, R)
    D, E = create_DE(lf, lf.SQ, lf.SQJ, R)
    new(lf, B, C, R, D, E)
  end
end

function create_A(lf::LurieFormAd)
  @import_fields lf.sys
  E = sparse([1:np;1:np], [from; to], [ones(np); -ones(np)], np, nb)

  A = [ spzeros(npq, npv+ npq)    speye(npq) #voltage pertubation limits
        spzeros(np, npq)   0*E[:,pv]  E[:,pq] #voltage discrepancy limits
        E[:,pv]  E[:,pq]   spzeros(np, npq) ] #angular discrepancy limits
  return A
end

# assuming δs = R*u
function create_BC(lf::LurieFormAd, A::AbstractMatrix, R::AbstractMatrix)
  B =  A*(lf.J⁻¹)*R
  C = -A*(lf.J⁻¹)*lf.S
  return (B, C)
end

function create_DE(lf::LurieFormAd, SQ::AbstractMatrix, SQJ::AbstractMatrix, R::AbstractMatrix)
  D = SQJ*(lf.J⁻¹)*R
  E = SQ - SQJ*(lf.J⁻¹)*lf.S
  return (D, E)
end

function create_box_R(lf::LurieFormAd, r::Real)
  @import_fields lf.sys

  # Find loads with nonzero value
  nzpq = find(b->(p⁰[pq[b]]!=0)||(q⁰[pq[b]]!=0), 1:npq)
  rnzpq = r*ones(length(nzpq))

  return create_loadability_R(lf, nzpq, rnzpq)
end

create_loadability_R(lf::LurieFormAd) = create_loadability_R(lf, lf.s⁰)
create_loadability_R(lf::LurieFormAd, ds::Vector) = reshape(ds, length(ds), 1)

function create_loadability_R(lf::LurieFormAd, nzpq::Vector, r::Vector)
  @import_fields lf.sys
  v⁰ = lf.sys.v⁰

  R = zeros(npv + 2*npq, length(nzpq))
  for b in 1:length(nzpq)
    R[npv + nzpq[b], b] = r[b]*p⁰[pq[nzpq[b]]]/(v⁰[pq[nzpq[b]]]^2)
    R[npv + npq + nzpq[b], b] = r[b]*q⁰[pq[nzpq[b]]]/(v⁰[pq[nzpq[b]]]^2)
  end

  return R
end

function create_RX_branch(cert::BaseLurieAdCertificate)
    sys = cert.lf.sys
    np = cert.lf.sys.np
    RX_branch = zeros(np,1)*im
    for p in 1:np
      RX_branch[p] = -1/sys.Y[sys.from[p],sys.to[p]]
    end
    return real(RX_branch), imag(RX_branch), RX_branch

end

function create_Imax(lf::LurieFormAd, rate::Real)
  YF = lf.sys.Yᶠ' - lf.sys.Yᵗ'
  V⁰ = lf.sys.v⁰ .* exp(lf.sys.ϑ⁰*im)
  Iₑ = YF*V⁰

  return abs(Iₑ)*rate
end

function get_Qlims(cert::BaseLurieAdCertificate, pm::GenericPowerModel)
  q⁰ = cert.lf.sys.q⁰
  net = cert.lf.sys.buses
  pv = cert.lf.sys.pv
  sl = cert.lf.sys.sl
  npv = length(pv)
  Qmin = []
  Qmax = []
  genbus = []
  for (gen, gen_data) in pm.data["gen"]
    Qmin = [Qmin; gen_data["qmin"]]
    Qmax = [Qmax; gen_data["qmax"]]
    genbus = [genbus; gen_data["gen_bus"]]
  end
   ngen = length(genbus)
   pos = []
  for k in 1:ngen
      pos = [pos; find(net -> net == genbus[k], net)]
  end
  pos[find(pos -> pos == sl[1], pos)] = 1e+7

  #reorder
  pos_in_pm = sortperm(pos)
  Qmin = Qmin[pos_in_pm]
  Qmax = Qmax[pos_in_pm]

  #allowed deviations
  dQmin = 1e+7*ones(npv)
  dQmax = 1e+7*ones(npv)

  for k in 1:npv
    if Qmin[k] <= q⁰[pv[k]]
       dQmin[k] = q⁰[pv[k]] - Qmin[k]
    end
    if q⁰[pv[k]]<= Qmax[k]
       dQmax[k] = -q⁰[pv[k]] + Qmax[k]
     end
  end
  return [dQmax dQmin]
end

function cal_box(R, Δ₁u, Δ₁v, npv, buspJ, basep, v⁰)
    resB = zeros(2, 2)
    res = R*Δ₁u[:,1]
    for k in 1:length(buspJ)
        resB[:, k] = -res[[npv+buspJ[1],npv+buspJ[2]]]
    end
    for k in 1:length(buspJ)
      # resB[k, 1] = exp(log(v⁰[buspJ[k]]) - Δ₁v[buspJ[k],2])^2 * (resB[k, 1] + basep[k]/(v⁰[buspJ[k]])^2) #max
      # resB[k, 2] = exp(log(v⁰[buspJ[k]]) - Δ₁v[buspJ[k],2])^2 * (-resB[k, 2] + basep[k]/(v⁰[buspJ[k]])^2) #min
      resB[k, 1] = (v⁰[buspJ[k]] - Δ₁v[buspJ[k],2])^2 * (resB[k, 1] + basep[k]/(v⁰[buspJ[k]])^2) #max
      resB[k, 2] = (v⁰[buspJ[k]] - Δ₁v[buspJ[k],2])^2 * (-resB[k, 2] + basep[k]/(v⁰[buspJ[k]])^2) #min
    end

    return resB
end

function polytope_plot(casename,busp,nump,resB,basep,nb,ratio,thermalrate,dV,stepsize,adaptive,limitQQ,limitIQ,Δ₁u)
    if Δ₁u[1,1] > 0
        @mput casename busp nump resB basep nb ratio thermalrate dV stepsize adaptive limitQQ limitIQ
        eval_string("
                    close all
                    addpath(genpath('/Users/hunghtd/Dropbox (Personal)/Journal/Feasibility region/matpower6.0b2'));

                    resR = constrainedCPF(casename, busp, nump, thermalrate, dV, stepsize, adaptive, limitQQ, limitIQ, 2*pi);

                    polytope_plot(resB, resR, basep, busp)

                    ")
    end

end

function add_Δ₁y_constraints_LP(cert::BaseLurieAdCertificate, m::Model, thermalrate::Real, T, K, dV::Real, limitIQ)
  @import_fields cert.lf.sys

  Δ₁y = getvariable(m, :Δ₁y)

  #current limits
  Imax = create_Imax(cert.lf, 2)
  R, X, Z = create_RX_branch(cert)
  Z = abs(Z)
  v⁰ = cert.lf.sys.v⁰
  ϑ⁰ = cert.lf.sys.ϑ⁰
  vt⁰ = v⁰[to]
  ρₑ⁰ = log(v⁰[from]./v⁰[to])
  Θₑ⁰ = ϑ⁰[from] - ϑ⁰[to]
  npvpq = npv+npq

  @variables m begin
    Δ₁v[1:npq,±] >= 0(start = 0.0)
    Δ₁ρ[1:np,±] >= 0(start = 0.0)
    Δ₁Θ[1:np,±] >= 0(start = 0.0)

    Δ₁ρᵐ[1:np] >= 0, (start = 0.0) #max of   Δ₁ρ\^±
    Δ₁Θᵐ[1:np] >= 0, (start = 0.0)
  end

  @constraints m begin
    [b in 1:npq,   s in ±], Δ₁y[b,s] == Δ₁v[b,s]
    [b in 1:npq,   s in ±], Δ₁y[b,s] <= dV
    # [p in 1:np,    s in ±], Δ₁y[npq + p,s] == Δ₁ρ[p,s]
    [p in 1:np,    s in ±], Δ₁y[npq + p,s] == Δ₁ρᵐ[p]
    # [p in 1:np,    s in ±], Δ₁y[npq + p,s] <= dV
    # [p in 1:np,    s in ±], Δ₁y[npq + np + p,s] == Δ₁Θ[p,s]
    [p in 1:np,    s in ±], Δ₁y[npq + np + p,s] == Δ₁Θᵐ[p]
    # [p in 1:np,    s in ±], Δ₁y[npq + np + p,s] <= T

    [b in 1:npq,   s in ±], Δ₁Θ[b,s] <= Δ₁Θᵐ[b]
    [p in 1:np           ], Δ₁Θᵐ[p] <= T

    [p in 1:np,   s in ±], Δ₁ρ[p,s] <= Δ₁ρᵐ[p]
    [p in 1:np          ], Δ₁ρᵐ[p] <= K
  end

  if limitIQ > 0 #if current limits are considered
      for p in 1:np, s = ±
          t = to[p]
          s̄ = negate(s)
          istobuspq = length(find(pq -> pq == t,pq)) > 0
          dVmaxload = Z[p]*Imax[p]/(vt⁰[p]+dV) #max voltage drop %
          dVmaxgen = Z[p]*Imax[p]/vt⁰[p]

          if istobuspq #to bus is a pq bus 1.1 for 39, 1.065 for 57
              @constraints m begin
                (1+2*Δ₁ρᵐ[p])*exp(2*ρₑ⁰[p]) - 2*exp(ρₑ⁰[p])*cos(T+Θₑ⁰[p])*(T*(1+Δ₁ρᵐ[p]+ρₑ⁰[p]) + Δ₁Θᵐ[p]*(1+K+ρₑ⁰[p]) - T*(1+K+ρₑ⁰[p])) - 1.2 <= dVmaxload^2
              end
          else
              @constraints m begin
                (1+2*Δ₁ρᵐ[p])*exp(2*ρₑ⁰[p]) - 2*exp(ρₑ⁰[p])*cos(T+Θₑ⁰[p])*(T*(1+Δ₁ρᵐ[p]+ρₑ⁰[p]) + Δ₁Θᵐ[p]*(1+K+ρₑ⁰[p]) - T*(1+K+ρₑ⁰[p])) - 1.2 <= dVmaxgen^2
              end
          end
        end
    end

end

function add_Δ₂f_constraints_LP(cert::BaseLurieAdCertificate, m::Model, Qlim::AbstractMatrix, T, K, limitQQ)
  @import_fields cert.lf.sys
  coshcos = range_coshcos(cert.lf)
  sinhcos = range_sinhcos(cert.lf)
  coshsin = range_coshsin(cert.lf)
  sinhsin = range_sinhsin(cert.lf)

  Δ₂f = getvariable(m, :Δ₂f);
  Δ₁ρ = getvariable(m, :Δ₁ρ);
  Δ₁Θ = getvariable(m, :Δ₁Θ);
  Δ₁ρᵐ = getvariable(m, :Δ₁ρᵐ);
  Δ₁Θᵐ = getvariable(m, :Δ₁Θᵐ);
  Δ₁u = getvariable(m, :Δ₁u);

  #constraints on reactive powers
  pv = cert.lf.sys.pv
  v⁰ = cert.lf.sys.v⁰
  npv = length(pv)

  SQ = cert.lf.SQ #last rows of the full matrix S
  SQ⁺, SQ⁻ = split_matrix(SQ)

  B = cert.B
  B⁺,B⁻ = split_matrix(B)

  D = cert.D
  D⁺, D⁻ = split_matrix(D)

  E = cert.E
  E⁺, E⁻ = split_matrix(E)

  nu = size(B,2)
  nf = size(cert.C,2)

  @variables m begin
   Δ₁cos[1:np,⁻] >= 0.0, (start = 0.0)
   Δ₂cos[1:np,⁻] >= 0.0, (start = 0.0)
   Δ₁sin[1:np,±] >= 0.0, (start = 0.0)
   Δ₂sin[1:np,±] >= 0.0, (start = 0.0)
   Δ₁cosh[1:np,⁺] >= 0.0, (start = 0.0)
   Δ₂cosh[1:np,⁺] >= 0.0, (start = 0.0)
   Δ₁sinh[1:np,±] >= 0.0, (start = 0.0)
   Δ₂sinh[1:np,±] >= 0.0, (start = 0.0)
  end

  for p in 1:np
    f = from[p]
    t = to[p]
    @constraints m begin
          Δ₁cos[p,⁻] == (1-cos(T))*Δ₁Θᵐ[p]/T
          Δ₂cos[p,⁻] == (1-cos(T))*Δ₁Θᵐ[p]/T

          Δ₁cosh[p,⁺] == (cosh(K)-1)*Δ₁ρᵐ[p]/K
          Δ₂cosh[p,⁺] == (cosh(K)-1)*Δ₁ρᵐ[p]/K

          Δ₂f[coshcos[p],⁺] == Δ₂cosh[p,⁺]

          Δ₂f[coshcos[p],⁻] >= (cosh(K)-1)*(1-cos(T))*Δ₁Θᵐ[p]/T + Δ₂cos[p,⁻]
          Δ₂f[coshcos[p],⁻] >= (cosh(K)-1)*(1-cos(T))*Δ₁ρᵐ[p]/K + Δ₂cos[p,⁻]
    end

    for s in ±
      s̄ = negate(s)
      @constraints m begin
            Δ₁sin[p,s] == Δ₁Θ[p,s]
            Δ₂sin[p,s] == (T-sin(T))*Δ₁Θ[p,s̄]/T
            # Δ₁sin[p,s] == Δ₁Θᵐ[p]
            # Δ₂sin[p,s] == (T-sin(T))*Δ₁Θᵐ[p]/T

            Δ₁sinh[p,s] == sinh(K)*Δ₁ρ[p,s]/K
            Δ₂sinh[p,s] == (sinh(K)-K)*Δ₁ρ[p,s]/K
            # Δ₁sinh[p,s] == sinh(K)*Δ₁ρᵐ[p]/K
            # Δ₂sinh[p,s] == (sinh(K)-K)*Δ₁ρᵐ[p]/K

            Δ₂f[coshsin[p],s] >= (cosh(K)-1)*Δ₁Θ[p,s] + Δ₂sin[p,⁺]
            # Δ₂f[coshsin[p],s] >= (cosh(K)-1)*Δ₁Θᵐ[p] + Δ₂sin[p,⁺]
            Δ₂f[coshsin[p],s] >= (cosh(K)-1)*Δ₁ρᵐ[p]*T/K + Δ₂sin[p,⁺]

            Δ₂f[sinhcos[p],s] >= sinh(K)*(1-cos(T))*Δ₁Θᵐ[p]/T + Δ₂sinh[p,s]
            Δ₂f[sinhcos[p],s] >= sinh(K)*(1-cos(T))*Δ₁ρ[p,s̄]/K + Δ₂sinh[p,s]
            # Δ₂f[sinhcos[p],s] >= sinh(K)*(1-cos(T))*Δ₁ρᵐ[p]/K + Δ₂sinh[p,s]

            Δ₂f[sinhsin[p],s] >= sinh(K)*Δ₁Θᵐ[p]
            Δ₂f[sinhsin[p],s] >= sinh(K)*Δ₁ρᵐ[p]*T/K
      end
    end

    if limitQQ > 0
        for k in 1:npv, s in ±
          s̄ = negate(s)
          @constraints m begin
            Qlim[k,s]/(v⁰[pv[k]]^2) >=
                sum(D⁺[k,j]*Δ₁u[j,s] + D⁻[k,j]*Δ₁u[j,s̄] for j=1:nu) +
                sum(E⁺[k,j]*Δ₂f[j,s] + E⁻[k,j]*Δ₂f[j,s̄] for j=1:nf)
          end
        end
    end

  end
end

function add_Δ₁y_constraints_woLP(cert::BaseLurieAdCertificate, m::Model, dV::Real, thermalrate::Real, limitIQ)# Δ₁vᴼ, Δ₁ρᴼ, Δ₁Θᴼ)
  @import_fields cert.lf.sys

  Δ₁y = getvariable(m, :Δ₁y)

  #current limits
  Imax = create_Imax(cert.lf, thermalrate)
  R, X, Z = create_RX_branch(cert)
  Z = abs(Z)
  v⁰ = cert.lf.sys.v⁰
  ϑ⁰ = cert.lf.sys.ϑ⁰
  vt⁰ = v⁰[to]
  ρₑ⁰ = log(v⁰[from]./v⁰[to])
  Θₑ⁰ = ϑ⁰[from] - ϑ⁰[to]
  npvpq = npv+npq

  @variables m begin
    Δ₁v[1:npq,±] >= 0(start = 0.0)
    Δ₁ρ[1:np,±] >= 0(start = 0.0)
    Δ₁Θ[1:np,±] >= 0(start = 0.0)
  end

  @constraints m begin
    [b in 1:npq,   s in ±], Δ₁y[b,s] == Δ₁v[b,s]
    [b in 1:npq,   s in ±], Δ₁y[b,s] <= dV
    [p in 1:np,    s in ±], Δ₁y[npq + p,s] == Δ₁ρ[p,s]
    [p in 1:np,    s in ±], Δ₁y[npq + p,s] <= dV
    [p in 1:np,    s in ±], Δ₁y[npq + np + p,s] == Δ₁Θ[p,s]
    [p in 1:np,    s in ±], Δ₁y[npq + np + p,s] <= pi/3
  end

  if limitIQ > 0
        for p in 1:np, s = ±
            t = to[p]
            s̄ = negate(s)
            # sgn = (-1)^s̄
            istobuspq = length(find(pq -> pq == t,pq)) > 0
            dVmaxload = Z[p]*Imax[p]/(vt⁰[p]+dV) #max voltage drop %
            dVmaxgen = Z[p]*Imax[p]/vt⁰[p]

            if istobuspq #to bus is a pq bus
                @NLconstraints m begin
                    sqrt(exp(2*Δ₁ρ[p,⁺]+2*ρₑ⁰[p]) - 2*cos(Δ₁Θ[p,s]+Θₑ⁰[p])*exp(Δ₁ρ[p,⁺]+ρₑ⁰[p]) + 1) <= dVmaxload
                end
            else
                @NLconstraints m begin
                    sqrt(exp(2*Δ₁ρ[p,⁺]+2*ρₑ⁰[p]) - 2*cos(Δ₁Θ[p,s]+Θₑ⁰[p])*exp(Δ₁ρ[p,⁺]+ρₑ⁰[p]) + 1) <= dVmaxgen
                end
            end
          end
   end

end

function add_Δ₂f_constraints_woLP(cert::BaseLurieAdCertificate, m::Model, Qlim::AbstractMatrix, limitQQ)
  @import_fields cert.lf.sys

  Δ₂f = getvariable(m, :Δ₂f);
  Δ₁ρ = getvariable(m, :Δ₁ρ);
  Δ₁Θ = getvariable(m, :Δ₁Θ);
  Δ₁u = getvariable(m, :Δ₁u);

  #constraints on reactive powers
  pv = cert.lf.sys.pv
  v⁰ = cert.lf.sys.v⁰
  npv = length(pv)

  SQ = cert.lf.SQ #last rows of the full matrix S
  SQ⁺, SQ⁻ = split_matrix(SQ)

  B = cert.B
  B⁺,B⁻ = split_matrix(B)

  D = cert.D
  D⁺, D⁻ = split_matrix(D)

  E = cert.E
  E⁺, E⁻ = split_matrix(E)

  nu = size(B,2)
  nf = size(cert.C,2)

  coshcos = range_coshcos(cert.lf)
  sinhcos = range_sinhcos(cert.lf)
  coshsin = range_coshsin(cert.lf)
  sinhsin = range_sinhsin(cert.lf)

  @variables m begin
   Δ₁cos[1:np,±] >= 0.0(start = 0.0)
   Δ₂cos[1:np,⁻] >= 0.0(start = 0.0)
   Δ₁sin[1:np,±] >= 0.0(start = 0.0)
   Δ₂sin[1:np,±] >= 0.0(start = 0.0)
   Δ₁cosh[1:np,±] >= 0.0(start = 0.0)
   Δ₂cosh[1:np,⁺] >= 0.0(start = 0.0)
   Δ₁sinh[1:np,±] >= 0.0(start = 0.0)
   Δ₂sinh[1:np,±] >= 0.0(start = 0.0)
  end

  for p in 1:np
    f = from[p]
    t = to[p]
    @constraints m begin
        Δ₁cos[p,⁺] == 0.
        Δ₁cosh[p,⁻] == 0.

        Δ₂f[coshcos[p],⁺] == Δ₂cosh[p,⁺]
        Δ₂f[coshcos[p],⁻] == Δ₁cosh[p,⁺]*Δ₁cos[p,⁻] + Δ₂cos[p,⁻]
    end

    for s in ±
      s̄ = negate(s)
      @NLconstraints m begin
            Δ₁cos[p,⁻] >= 1. - cos(Δ₁Θ[p,s])
            Δ₂cos[p,⁻] >= 1. - cos(Δ₁Θ[p,s])

            Δ₁cosh[p,⁺] >= cosh(Δ₁ρ[p,s]) - 1.
            Δ₂cosh[p,⁺] >= cosh(Δ₁ρ[p,s]) - 1.
            Δ₂sin[p,s] == Δ₁Θ[p,s̄] - sin(Δ₁Θ[p,s̄])
            Δ₂sinh[p,s] == sinh(Δ₁ρ[p,s̄]) - Δ₁ρ[p,s̄]

            Δ₂f[coshsin[p],s] == Δ₁cosh[p,⁺]*sin(Δ₁Θ[p,s]) + Δ₁Θ[p,s̄] - sin(Δ₁Θ[p,s̄])
            Δ₂f[sinhcos[p],s] == sinh(Δ₁ρ[p,s̄])*Δ₁cos[p,⁻] + sinh(Δ₁ρ[p,s]) - Δ₁ρ[p,s]

            Δ₂f[sinhsin[p],s] >= sinh(Δ₁ρ[p,⁺])*sin(Δ₁Θ[p,s])
            Δ₂f[sinhsin[p],s] >= sinh(Δ₁ρ[p,⁻])*sin(Δ₁Θ[p,s̄])
      end
    end

   if limitQQ > 0
        for k in 1:npv, s in ±
          s̄ = negate(s)
          @constraints m begin
            Qlim[k,s]/(v⁰[pv[k]]^2) >=
                sum(D⁺[k,j]*Δ₁u[j,s] + D⁻[k,j]*Δ₁u[j,s̄] for j=1:nu) +
                sum(E⁺[k,j]*Δ₂f[j,s] + E⁻[k,j]*Δ₂f[j,s̄] for j=1:nf)
          end
        end
    end

  end
end

export BaseLurieAdCertificate, create_A, create_BC
export create_box_R, create_loadability_R
export create_RX_branch, create_Imax, get_Qlims, cal_box, polytope_plot
export add_Δ₁y_constraints_woLP, add_Δ₂f_constraints_woLP
