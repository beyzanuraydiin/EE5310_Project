include("auto_Ymatrix.jl")
include("power_flow_beyz.jl")

using Plots, LinearAlgebra

P_net = initialize_parameters()
Y_bus_E =  Ybus_By_IncidenceMatrix()

# Construct the abs(Ybus matrix)
P_net = P_net[setdiff(1:end, 21), :]
Y_bus_E = Y_bus_E[setdiff(1:end, 21), setdiff(1:end, 21)]

a = -imag(Y_bus_E)\P_net

plot((a))