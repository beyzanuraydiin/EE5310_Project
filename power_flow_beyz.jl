# %% 
using LinearAlgebra
include("PF_data.jl")
include("auto_Ymatrix.jl")
using Plots

# -------- initialize parameters ---------
# slack bus voltage and angle
function initialize_parameters()
    v0 = 1 #pu
    delta0 = 0 
    S_base = 100 #100 MVA
    # Load data
    nbus = length(bus_data[:,1])#number of buses
    bus_no    = bus_data[:,1]
    bus_type  = bus_data[:,2] #type(1=PQ,2=PV,3=SLACK)
    pv_buses  = findall(x -> x ==2, bus_type) #we know p and v - find delta and q
    pq_buses  = findall(x -> x ==1, bus_type) #we know p and q - find delta and v
    slack_bus = findall(x-> x==3, bus_type)
    no_pv = length(pv_buses) #number of pv buses
    no_pq = length(pq_buses)
    # bus_V   = 
    bus_Pload = zeros(nbus,1)
    bus_Qload = zeros(nbus,1) #create a column of zeros and store the P loads -*and Q loads
    
    for b in 1:length(load_data[:,1])
        bus_Pload[Int(load_data[b,1]),1] = load_data[b,3]/S_base      
    end
    
    for b in 1:length(load_data[:,1])
        if Int(load_data[b,1]) in pq_buses
            bus_Qload[Int(load_data[b,1]),1] = (load_data[b,3]*tan(acos.(load_data[b,2])))/S_base 
        end
    end
    
    bus_Pgi = zeros(nbus,1)
    bus_Qgi = zeros(nbus,1)
    # the sum of all generation in a bus
    for g in 1:nbus
        bus_Pgi[g,1] = sum(gen_data[gen_data[:, 1] .== g, 5])/S_base 
    end
    net_wind = zeros(nbus,1)
    for w in 1:length(wind_data[:,1])
        net_wind[Int(wind_data[w,1]),1] = wind_data[w,5]/S_base  
    end
    P_net = bus_Pgi + net_wind - bus_Pload
    Q_net = bus_Qgi - bus_Qload
    
    #remove the slack bus
    # P_net = P_net[setdiff(1:end, slack_bus), :]
    # Q_net = Q_net[setdiff(1:end, slack_bus), :]
    y = [P_net;Q_net]
    
    Y_matrix  = construct_admittance_matrix()
    
    
    G = real(Y_matrix)
    B = imag(Y_matrix)
    V = ones(nbus,1)
    for v in 1:nbus
        if !isempty(gen_data[gen_data[:, 1] .== v, 4])
            v_pv = gen_data[gen_data[:, 1] .== v, 1][1]
            if v_pv in pv_buses
                V[v,1] = gen_data[gen_data[:, 1] .== v, 4][1]
            end
        end
    end
    # remove slack bus
    # V = V[setdiff(1:end, slack_bus), :]
# ------------------- compute y and x matrices ---------------------
    #remove the slack bus from Y_matrix
    # Y_matrix = Y_matrix[setdiff(1:end, slack_bus), setdiff(1:end, slack_bus)]
    thetas = angle.(Y_matrix)

    d_x = ones(nbus,1)*10e-8
    d_x[21,1] = 0
    v_x = V
    v_x[21,1] = 1
    x_dv = [d_x; v_x]

    return P_net, Q_net, Y_matrix, V, d_x, v_x, nbus, pv_buses, pq_buses, bus_type
end


#----------------------------------------------------------------------
#compute the Jacobian matrix
# %%
function compute_jacobian(Y_matrix, nbus, v_x,d_x)
    P_net, Q_net, Y_matrix, V, d_x, v_x, nbus, pv_buses, pq_buses, bus_type = initialize_parameters()

    J1 = zeros(Complex, nbus,nbus)
    J2 = zeros(Complex, nbus,nbus)
    J3 = zeros(Complex, nbus,nbus)
    J4 = zeros(Complex, nbus,nbus)
    # off-diagonal entries
    for k in 1:nbus
        for j in 1:nbus
            if k !=j
                J1[k,j] = v_x[k]*abs(Y_matrix[k,j])*v_x[j]*sin(d_x[k]-d_x[j]-thetas[k,j])
                J2[k,j] = v_x[k]*abs(Y_matrix[k,j])*cos(d_x[k]-d_x[j]-thetas[k,j])  
                J3[k,j] = -v_x[k]*abs(Y_matrix[k,j])*v_x[j]*cos(d_x[k]-d_x[j]-thetas[k,j])
                J4[k,j] = v_x[k]*abs(Y_matrix[k,j])*sin(d_x[k]-d_x[j]-thetas[k,j])
            end
        end
    end

    y1_tot = zeros(Complex,nbus,1)
    y2_tot = zeros(Complex,nbus,1)
    # diagonal entries
    
    for k in 1:nbus
            for j in 1:nbus 
                if j != k
                    y1_tot[k] = y1_tot[k] + abs(Y_matrix[k,j])*v_x[j]*sin(d_x[k]-d_x[j]-thetas[k,j])
                end
            end
            for j in 1:nbus
                y2_tot[k] = y2_tot[k] + abs(Y_matrix[k,j])*v_x[j]*cos(d_x[k]-d_x[j]-thetas[k,j]) 
            end
            J1[k,k] = -v_x[k]*y1_tot[k]
            J2[k,k] = v_x[k]*abs(Y_matrix[k,k])*cos(thetas[k,k]) + y2_tot[k]
            J3[k,k] = v_x[k]*y1_tot[k]
            J4[k,k] = -v_x[k]*abs(Y_matrix[k,k])*sin(thetas[k,k]) + y1_tot[k]
        
    end
    Jac = [J1 J2;J3 J4]

    return Jac
end
# %%
function NewtonRaphson(d_x, v_x)
    P_net, Q_net, Y_matrix, V, d_x, v_x, nbus, pv_buses, pq_buses, bus_type = initialize_parameters()

    b = 0
    while true
        F = compute_P_Q(v_x, d_x,Y_matrix,thetas,P_net,Q_net)
        b = b + 1
        Jac = compute_jacobian(Y_matrix, nbus, v_x,d_x)

        if det((Jac)) < eps() #Machine Epsilon for Float64
            println("Jacobian is non-invertible.")
            break
        end
        X_ans = [d_x; v_x]
        X_ans = X_ans - Jac \ F
        
        d_x = d_x + X_ans[1:nbus,1]
        v_x = v_x + X_ans[nbus+1:end,1]
        
        if b ==5
            break
        end
    end 

end
# %%
function compute_P_Q(v_x,d_x,Y_matrix,thetas,P_net,Q_net)
    P_net, Q_net, Y_matrix, V, d_x, v_x, nbus, pv_buses, pq_buses, bus_type = initialize_parameters()

    P = zeros(Complex,nbus,1)
    Q = zeros(Complex,nbus,1)
    for k in 1:nbus
        for n in 1:nbus
            P[k] = P[k] + v_x[k]*(abs(Y_matrix[k,n])*v_x[n]*cos(d_x[k]-d_x[n]-thetas[k,n])) #y_k
            Q[k] = Q[k] + v_x[k]*(abs(Y_matrix[k,n])*v_x[n]*sin(d_x[k]-d_x[n]-thetas[k,n])) #y_k+N   
        end
    end 
    Px = P - P_net
    Qx = Q - Q_net
    F = [Px; Qx]
    return F
end
# %%
