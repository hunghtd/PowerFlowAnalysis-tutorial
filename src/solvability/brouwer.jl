using Gurobi
using HSL
⁺ = 1; ⁻ = 2; ± = 1:2
negate(i::Integer) = i==⁺ ? ⁻ : ⁺
split_matrix(C) = (max.(C, 0.0), max(-C, 0.0))

# Describes the system Δ₁y = AΔ₁x = (B⁺ - B⁻)Δ₁u + (C⁺ - C⁻)Δ₂{f}
function brouwer_certificate_LP(sys, B, C, Qlim, thermalrate, T, K, dV, Iter, limitQQ, limitIQ;
  region = :default, objective = :default)
  ny = size(B,1)
  nu = size(B,2)
  nf = size(C,2)

  B⁺,B⁻ = split_matrix(B)
  C⁺,C⁻ = split_matrix(C)

  m = Model(solver = GurobiSolver(Method=2,Crossover=0, NodefileStart=0.1))#,outlev=1))

  @variables m begin
    Δ₁y[1:ny,±] >= 0.0, (start = 0.0)
    Δ₂f[1:nf,±] >= 0.0, (start = 0.0)
    Δ₁u[1:nu,±] >= 0.0, (start = 0.0)
  end
  for k in 1:ny, s in ±
    s̄ = negate(s)
    @constraints m begin
      Δ₁y[k,s] >=
        sum(B⁺[k,j]*Δ₁u[j,s] + B⁻[k,j]*Δ₁u[j,s̄] for j=1:nu) +
        sum(C⁺[k,j]*Δ₂f[j,s] + C⁻[k,j]*Δ₂f[j,s̄] for j=1:nf)
    end
  end
  add_Δ₁y_constraints_LP(sys, m, thermalrate, T, K, dV, limitIQ)
  add_Δ₂f_constraints_LP(sys, m, Qlim, T, K, limitQQ)

  if region == :cone
    @constraint(m, [k in 1:nu], Δ₁u[k,⁻] == 0.0)
  elseif region == :uniform
    @constraint(m, [k in 1:nu, s in ±], Δ₁u[k,s] == Δ₁u[1,⁺])
  elseif region == :custom
    add_Δ₁u_constraints(m, sys)
  end # otherwise, don't impose any additional constraints

  if objective == :custom
    add_objective(m, sys)
  elseif objective == :gauss
    s2 = √2
    @NLexpression(m, lpr[k in 1:nu], log(erf(Δ₁u[k,⁺]/s2) + erf(Δ₁u[k,⁻]/s2)))
    @variable(m, logprob)
    @NLconstraint(m, logprob == sum(lpr[k] for k in 1:nu) - 0*nu*log(2))
    @objective(m, Max, logprob)
  elseif objective == :default
    @objective(m, Max, Δ₁u[1,⁺])
  end

  return m
end

export brouwer_certificate_LP
