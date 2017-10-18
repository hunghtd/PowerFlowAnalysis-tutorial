type Network
  pairs::Vector{Tuple{Int,Int}}
  pair_idx::Dict{Tuple{Int,Int}, Int}
  np::Int

  buses::Vector{Int}
  bus_idx::Dict{Int, Int}
  nb::Int

  from::Vector{Int}
  to::Vector{Int}

  function Network(pairs::Vector{Tuple{Int,Int}}, buses::Vector{Int})
    pair_idx = Dict(p=>j for (j,p) in enumerate(pairs))
    np = length(pairs)

    bus_idx = Dict(b=>j for (j,b) in enumerate(buses))
    nb = length(buses)

    from  = [bus_idx[p[1]] for p in pairs]
    to  =   [bus_idx[p[2]] for p in pairs]

    new(pairs, pair_idx, np, buses, bus_idx, nb, from, to)
  end
end

function Network(pm::GenericPowerModel)
  pairs = [k for k in keys(pm.ref[:buspairs])]
  buses = [k for k in keys(pm.ref[:bus])]

  Network(pairs, buses)
end

function calc_YYY(pm::GenericPowerModel, net::Network)
  Yᵈ = zeros(Complex, net.nb);
  Yᶠ = zeros(Complex, net.np);
  Yᵗ = zeros(Complex, net.np);

  for (br, br_data) in pm.ref[:branch]
    g, b = PowerModels.calc_branch_y(br_data)
    tr, ti = PowerModels.calc_branch_t(br_data)

    Ys = br_data["br_status"]*(g + im*b)
    Ytt = Ys + im*br_data["br_status"]*br_data["br_b"]/2

    br_ft = (br_data["f_bus"], br_data["t_bus"])
    f_idx = net.bus_idx[br_ft[1]]
    t_idx = net.bus_idx[br_ft[2]]

    Yᵈ[f_idx] = Yᵈ[f_idx] + Ytt/br_data["tap"]^2
    Yᵈ[t_idx] = Yᵈ[t_idx] + Ytt

    if haskey(net.pair_idx, br_ft)
      br_idx = net.pair_idx[br_ft]
      Yᶠ[br_idx] = Yᶠ[br_idx] - Ys/(tr - im*ti)
      Yᵗ[br_idx] = Yᵗ[br_idx] - Ys/(tr + im*ti)
    elseif haskey(net.pair_idx, reverse(br_ft))
      br_idx = net.pair_idx[reverse(br_ft)]
      Yᵗ[br_idx] = Yᵗ[br_idx] - Ys/(tr - im*ti)
      Yᶠ[br_idx] = Yᶠ[br_idx] - Ys/(tr + im*ti)
    else
      ErrorException("Inconsistent pair and branch data")
    end
  end

  for (bus, bus_data) in pm.ref[:bus]
    Yᵈ[net.bus_idx[bus]] = Yᵈ[net.bus_idx[bus]] + bus_data["gs"] + im*bus_data["bs"]
  end

  Yᵈ = spdiagm(Yᵈ, (0))
  Yᶠ = sparse(net.from, 1:net.np, Yᶠ, net.nb, net.np)
  Yᵗ = sparse(net.to,   1:net.np, Yᵗ, net.nb, net.np)

  return (Yᵈ, Yᶠ, Yᵗ)
end

function get_bus_types(pm::GenericPowerModel, net::Network)
  pq    = find(b -> pm.ref[:bus][b]["bus_type"] == 1, net.buses)
  pv    = find(b -> pm.ref[:bus][b]["bus_type"] == 2, net.buses)
  slack = find(b -> pm.ref[:bus][b]["bus_type"] == 3, net.buses)

  return (pq, pv, slack)
end

function get_pq(pm::GenericPowerModel, net::Network)
  p = [-pm.ref[:bus][idx]["pd"]  for idx in net.buses]
  q = [-pm.ref[:bus][idx]["qd"]  for idx in net.buses]

  for data in values(pm.ref[:gen])
    idx = net.bus_idx[data["gen_bus"]]
    p[idx] = p[idx] + data["gen_status"]*data["pg"]
    q[idx] = q[idx] + data["gen_status"]*data["qg"]
  end
  return (p, q)
end

function get_vϑ(pm::GenericPowerModel, net::Network)
  v = [pm.ref[:bus][idx]["vm"]        for idx in net.buses]
  ϑ = [pm.ref[:bus][idx]["va"]*pi/180 for idx in net.buses]
  return (v, ϑ)
end


function get_sol_vϑ(sol::Dict, net::Network)
  v = [sol["solution"]["bus"][string(idx)]["vm"]        for idx in net.buses]
  ϑ = [sol["solution"]["bus"][string(idx)]["va"]*pi/180 for idx in net.buses]
  return (v, ϑ)
end

# function zero_loading(pm::GenericPowerModel)
#     pm0 = pm;
#     allbuses = keys(pm.data["bus"])
#     allgens = keys(pm.data["gen"])
#
#     for idx in allbuses
#       pm0.data["bus"][idx]["pd"] = 0.01*pm.data["bus"][idx]["pd"]
#       pm0.data["bus"][idx]["qd"] = 0.01*pm.data["bus"][idx]["qd"]
#     end
#     for idx in allgens
#       pm0.data["gen"][idx]["pg"] = 0.01*pm.data["gen"][idx]["pg"]
#     end
#    return pm0
# end

function zero_loading(case::Dict)
    case0 = case;
    allbuses = keys(case["bus"])
    allgens = keys(case["gen"])

    for idx in allbuses
      case0["bus"][idx]["pd"] = 0.01*case["bus"][idx]["pd"]
      case0["bus"][idx]["qd"] = 0.01*case["bus"][idx]["qd"]
    end
    for idx in allgens
      case0["gen"][idx]["pg"] = 0.01*case["gen"][idx]["pg"]
    end
   return case0
end

export Network, calc_YYY, get_bus_types, get_pq, get_vϑ, get_sol_vϑ, zero_loading
