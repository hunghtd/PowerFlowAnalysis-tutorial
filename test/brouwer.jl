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
# ratio = pqload[buspJ[1],1]/pqload[buspJ[2],1]
ϕ = acot(ratio)
resB = zeros(2, 2)
v⁰ = lf.sys.v⁰
R = create_loadability_R(lf, buspJ, [cos(ϕ), sin(ϕ)])

lc = BaseLurieAdCertificate(lf, R)

Qlim = get_Qlims(lc, pm)

dV = 0.05
thermalrate = 2.

mNP = brouwer_certificate_woLP(lc, lc.B, lc.C, 200, Qlim, dV, thermalrate; region=:uniform)
# DegeneracyHunter.printVariableDiagnostics(mNP)
sol = solve(mNP)

DegeneracyHunter.printVariableDiagnostics(mNP)
Δ₁u = getvalue(getvariable(mNP, :Δ₁u))
Δ₁v = getvalue(getvariable(mNP, :Δ₁v))

#options for CPF
# dV = findmax(Δ₁v[1,:])[1]
# dV = 0.1
stepsize = 0.05
adaptive = 0

basep = -sys.p⁰[sys.pq[buspJ]]

res = R*Δ₁u[:,1]
for k in 1:2
    resB[:, k] = -res[[npv+buspJ[1],npv+buspJ[2]]]
end
for k in 1:length(buspJ)
  resB[k, 1] = exp(log(v⁰[buspJ[k]]) - Δ₁v[buspJ[k],2])^2 * (resB[k, 1] + basep[k]/(v⁰[buspJ[k]])^2) #max
  resB[k, 2] = exp(log(v⁰[buspJ[k]]) - Δ₁v[buspJ[k],2])^2 * (-resB[k, 2] + basep[k]/(v⁰[buspJ[k]])^2) #min
end

# resB =  -0.0788118
#  -0.044624

nump = 4
@mput casename busp nump resB basep nb ratio thermalrate dV stepsize adaptive
eval_string("close all;
            addpath(genpath('/Users/hunghtd/Documents/MATLAB/matpower6.0b2'));
            %resR = CPFfun_Imax3(casename, busp, nump, 2*pi);
            %resR = CPFfun_IQmax(casename, busp, nump, dV, 2*pi);
            %resR = CPFfun_Imax_noQ(casename, busp, nump, dV, 2*pi);
            resR = CPFfun_IQmax_new(casename, busp, nump, thermalrate, dV, stepsize, adaptive, 2*pi);

            A = [1 0; 0 1; -1 0; 0 -1];

            b = [resB(1, 1); resB(2, 1); -resB(1, 2); -resB(2, 2)];
            V = con2vert(A, b);
            x = V(:,1);
            y = V(:,2);
            h = convhull(x,y);

            plot(x(h)+0*basep(1), y(h)+0*basep(2), 'linewidth', 3)

            hold on
            plot(resR(:,1), resR(:,2), 'r','linewidth', 2.5)
            scatter(basep(1), basep(2), [],'k', '*', 'LineWidth',4)
            %axis([basep(1) inf basep(2) inf])
            %title(['Solvability boundaries for IEEE ',num2str(nb),' bus system'])
            xlabel(['Active power at load bus ', num2str(busp(1)), ' (P_{', num2str(busp(1)), '})'])
            ylabel(['Active power at load bus ', num2str(busp(2)), ' (P_{', num2str(busp(2)), '})'])
            set(gca,'Fontname','Times New Roman','fontsize',24)
            legend( 'Inner approximation', 'Feasibility boundary', 'Base operating point')
              "
              )
