include("../src/PowerFlowAnalysis.jl")

using JuMP
using Ipopt
using PowerModels
using PowerFlowAnalysis
using MATLAB
using PyPlot

casename = "case9" #Matpower test case
case = PowerModels.parse_file(join(["data/", casename, ".m"]))
pm = ACPPowerModel(case)
sol = run_ac_pf(case, IpoptSolver())
sys = PowerFlowSystem(pm, sol)
lf = LurieFormAd(sys)

# Choose between different choices of the u - constraints
net = Network(pm)
pqload = sys.p⁰[sys.pq]
nzload = find(pqload -> pqload < 0, pqload)
buspJ = nzload[[1,2]]##the bus ordering in this code
busp = [net.buses[k] for k in sys.pq[buspJ]] #the bus ordering in Matpower cases

npv = length(sys.pv)
nb = sys.nb

ratio = 1; #aspect ratio of two loads 
ϕ = acot(ratio)
v⁰ = lf.sys.v⁰
R = create_loadability_R(lf, buspJ, [cos(ϕ), sin(ϕ)]) #loading scenario

lc = BaseLurieAdCertificate(lf, R)

Qlim = get_Qlims(lc, pm)
dV = 0.01 #voltage derivation limit
thermalrate = 2.
T = 0.01#angle separation limit
K = 0.01#logarithmic voltage difference limit

limitQQ = 1 #if reactive power constraints are considered
limitIQ = 1 #if reactive thermal constraints are considered

#core optimization problem for construction
mLP = brouwer_certificate_LP(lc, lc.B, lc.C, Qlim, thermalrate, T, K, dV, 100, limitQQ, limitIQ; region=:uniform)
sol_LP = solve(mLP)

Δ₁u = getvalue(getvariable(mLP, :Δ₁u))
Δ₁v = getvalue(getvariable(mLP, :Δ₁v))

#options for Continuation Power Flow
dV = findmax(Δ₁v[1,:])[1]
stepsize = 0.03
adaptive = 0

basep = -sys.p⁰[sys.pq[buspJ]] #the base point
resB = cal_box(R, Δ₁u, Δ₁v, npv, buspJ, basep, v⁰) #inner approx. region
nump = 4 #number of points to plot

polytope_plot(casename,busp,nump,resB,basep,nb,ratio,thermalrate,dV,stepsize,adaptive,limitQQ,limitIQ,Δ₁u)
