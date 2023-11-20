# %% 
using LinearAlgebra
include("PF_data.jl")
using Plots

# STEP 0
# %%
# Extract data from the IEEE RTS 96 FullSystemData
function initialize_data()
    n = length(bus_data[:,1]) #number of nodes/buses of the system
    frombus = Int.(real(line_data[:,1])) 
    tobus   = Int.(real(line_data[:,2]))#
    l = length(frombus)
    
    R = real(line_data[:,3]) 
    X = real(line_data[:,4])
    Z = R + X*im #impedance calculation using the given data
    
    ONTR = line_data[:,6]
    Gsh  = bus_data[:,3]     #conductance of shunt adm (Gsh,pu)
    Bsh  = bus_data[:,4]     #susceptance of shunt adm (Bsh,pu)
    B_shunt_ij = line_data[:,5] 

    return n,frombus,tobus,l,R,X,Z,ONTR,Gsh,Bsh,B_shunt_ij
end
function construct_admittance_matrix()
# %%
n,frombus,tobus,l,R,X,Z,ONTR,Gsh,Bsh,B_shunt_ij = initialize_data()
# %%
# STEP 1
Y_matrix  = zeros(ComplexF64,n,n);   #off-diagonal elements of Y bus matrix
Y_diag    = zeros(ComplexF64,n,1);     #diagonal elements of Y bus matrix
slack_bus = bus_data[bus_data[:, 2] .== 3, 1][1]

# %% ------------------------------------------
#time to add tap ratios and phase shifters
for i in 1:l
    if ONTR[i] == 0
        c = 1
    else
        c = ONTR[i]
    end
    y_t_ij = 1/Z[i]
    Y_matrix[(frombus[i]), (tobus[i])] += y_t_ij/conj(c)  #Yij = -y_ij/a*ij
    Y_matrix[(tobus[i]),(frombus[i])] += y_t_ij/c         #Yji = -y_ij/a_ij
end

for i in 1:l
    if ONTR[i] == 0
        c = 1
    else
        c = ONTR[i]
    end
    y_t_ij = 1/Z[i]
    Y_diag[(frombus[i])] = Y_diag[(frombus[i])] + y_t_ij./(abs.(c)^2) #Yii = Yii + y_ij/a^2
    Y_diag[(tobus[i])] = Y_diag[(tobus[i])] + y_t_ij                #Yjj = Yjj + y_ij
end


#time to add shunts -  they only affect diagonal elements
for i in 1:l
    shunt_tot = sum(line_data[line_data[:, 1] .== i, 5]) #sum of shunts connected to bus i
    
    Y_diag[(frombus[i])] = -Y_diag[(frombus[i])] +Gsh[(frombus[i])] +im*Bsh[(frombus[i])] + im*shunt_tot/2  #Yii
    Y_diag[(tobus[i])] = -Y_diag[(tobus[i])]  +Gsh[(frombus[i])] +im*Bsh[(frombus[i])]+ im*shunt_tot/2                #Yjj

   end

for i in 1:n
    Y_matrix[i,i] = Y_diag[i];
end

return Y_matrix
# %%
end
# %% -------------- C -----------------
Y_bus =  construct_admittance_matrix()
spy(abs.(Y_bus))
# %% -------------- D -----------------# print last 5 rows and columns
last_five = Y_bus[69:73, 69:73]
# %% -------------- E -----------------
#set all tap ratios to 1 -phase shifters to 0 
function Ybus_By_IncidenceMatrix()
    n,frombus,tobus,l,R,X,Z,ONTR,Gsh,Bsh,B_shunt_ij = initialize_data()

    l_d = ones(length(line_data[:,6]),1)
    # construct the incidence matrix  
    m = length(line_data[:,1]) #number of lines
    n = length(bus_data[:,1])
    E = zeros(m,n) #number of branches x number of buses
    for i in 1:m
        for j in 1:n
            if frombus[i] == j 
                E[i, j] = 1
            elseif tobus[i] == j
                E[i, j] = -1
            end
        end
    end
    Y_l = diagm(1 ./Z)
    Y_s = zeros(Complex,n,n)  
    for i in 1:m
        for j in 1:n
            if frombus[i] == j
                Y_s[j,j] = Y_s[j,j] + 0.5*im*B_shunt_ij[i]
            elseif tobus[i] == j
                Y_s[j,j] = Y_s[j,j] + 0.5*im*B_shunt_ij[i]
            end
        end
    end
    Y_bus_E = E' * Y_l * E + Y_s
    return Y_bus_E
end 
# %% -------------- F -----------------
Y_bus_E =  Ybus_By_IncidenceMatrix()
# spy(abs.(Y_bus_E))
Ybus_norm =  norm(Y_bus)
YbusE_norm = norm(Y_bus_E)
difff = Ybus_norm - YbusE_norm
if difff==0
    println("Norms are equal!")
else
    println("Norms are not equal!")
end

