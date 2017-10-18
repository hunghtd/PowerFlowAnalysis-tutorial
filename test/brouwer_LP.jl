include("../src/PowerFlowAnalysis.jl")
include("../src/degenhunter/DegeneracyHunter.jl")

using JuMP
using Ipopt
using PowerModels
using PowerFlowAnalysis
using MATLAB
using PyPlot

casename = "case300"
case = PowerModels.parse_file(join(["data/", casename, ".m"]))
pm = ACPPowerModel(case)
sol = run_ac_pf(case, IpoptSolver())
sys = PowerFlowSystem(pm, sol)
lf = LurieFormAd(sys)

# Choose between different choices of the u - constraints
net = Network(pm)
pqload = sys.p⁰[sys.pq]
nzload = find(pqload -> pqload < 0, pqload)
buspJ = nzload[[2,3]]#2, 3 for 300
busp = [net.buses[k] for k in sys.pq[buspJ]]

npv = length(sys.pv)
nb = sys.nb

ratio = 1;
# ratio = pqload[buspJ[2],1]/pqload[buspJ[1],1]
ϕ = acot(ratio)
v⁰ = lf.sys.v⁰
R = create_loadability_R(lf, buspJ, [cos(ϕ), sin(ϕ)])

lc = BaseLurieAdCertificate(lf, R)

Qlim = get_Qlims(lc, pm)
dV = 0.01 #0.002 for 39; 0.05 for 57; 0.01 for 300
thermalrate = 2.
T = 0.01#0.001 for 39, 0.015 for 118; 0.015 for pnnl300
K = 0.01#0.001 for 118; 0.1 for pnnl300

limitQQ = 1 #if reactive power constraints are considered
limitIQ = 1 #if reactive thermal constraints are considered

mLP = brouwer_certificate_LP(lc, lc.B, lc.C, Qlim, thermalrate, T, K, dV, 100, limitQQ, limitIQ; region=:uniform)
sol_LP = solve(mLP)

DegeneracyHunter.printVariableDiagnostics(mLP)
Δ₁u = getvalue(getvariable(mLP, :Δ₁u))
Δ₁v = getvalue(getvariable(mLP, :Δ₁v))

#options for CPF
dV = findmax(Δ₁v[1,:])[1]
stepsize = 0.03
adaptive = 0

basep = -sys.p⁰[sys.pq[buspJ]]
resB = cal_box(R, Δ₁u, Δ₁v, npv, buspJ, basep, v⁰)
nump = 4

polytop_plot(casename,busp,nump,resB,basep,nb,ratio,thermalrate,dV,stepsize,adaptive,limitQQ,limitIQ,Δ₁u)

# if Δ₁u[1,1] > 0
#     # polytop_plot_new(casename,busp,nump,resB,basep,nb,ratio,thermalrate,dV,stepsize,adaptive,limitQQ,limitIQ, Δ₁u)
#     @mput casename busp nump resB basep nb ratio thermalrate dV stepsize adaptive limitQQ limitIQ
#     eval_string("
#                 close all
#                 addpath(genpath('/Users/hunghtd/Dropbox (Personal)/Journal/Feasibility region/matpower6.0b2'));
#
#                 resR = constrainedCPF(casename, busp, nump, thermalrate, dV, stepsize, adaptive, limitQQ, limitIQ, 2*pi);
#
#                 polytop_plot(resB, resR, basep, busp)
#
#                 ")
# end
#
# resB
# @mget(resR)
