%% Newton Rapson Method Load Flow
% Beyzanur Aydin
clear all; clc;
%%Ybus construction
frombus = [1; 1; 2];
tobus = [2; 3; 3];
R = [0.05; 0.05; 0.05];
X = [0.1; 0.1; 0.1];
bus_no = [1; 2; 3];
Pd = [100; -180; 80];
Sbase = 100; % 100 MVA base
Pd_u = Pd / Sbase;
Qd = [0; 0; 0];
Qd_u = Qd / Sbase;
Z = R + 1i * X;
Y = 1 ./ Z;

nbus = length(bus_no); % number of buses
nline = length(frombus); % number of branches

% Construct the abs(Ybus matrix)
Y_od = zeros(nbus, nbus); % off-diagonal elements
for i = 1:nline
    Y_od(frombus(i), tobus(i)) = -Y(i);
    Y_od(tobus(i), frombus(i)) = -Y(i);
end

Y_diag = zeros(nbus, 1);
for i = 1:nbus
    Y_diag(frombus(i)) = Y_diag(frombus(i)) + 1 / Z(i);
    Y_diag(tobus(i)) = Y_diag(tobus(i)) + 1 / Z(i);
end

Ybus = diag(Y_diag) + Y_od;
Gsh = real(Ybus);
Bsh = imag(Ybus);
thetas = angle(Ybus);

% Defining the known variables
d1 = 0;
v1 = 1.0;
v3 = 1.0;

% Defining the power injections
p2 = Pd_u(2,1);
q2 = Qd_u(3,1);
p3 = Pd_u(3,1);

% Defining initial values of unknowns
d2 = 0.0;
d3 = 0.0;
v2 = 1.0;

X = [d2; d3; v2];

y=abs(Ybus);

syms d2 d3 v2
P2=y(2,1)*v1*v2*cos(thetas(2,1)+d1-d2)+y(2,2)*v2^2*cos(thetas(2,2))+y(2,3)*v3*v2*cos(thetas(2,3)+d3-d2)-p2;
P3=y(3,1)*v1*v3*cos(thetas(3,1)+d1-d3)+y(3,2)*v3*v2*cos(thetas(3,2)+d2-d3)+y(3,3)*v3^2*cos(thetas(3,3))-p3;
Q2=-y(2,1)*v1*v2*sin(thetas(2,1)+d1-d2)-y(2,2)*v2^2*sin(thetas(2,2))-y(2,3)*v3*v2*sin(thetas(2,3)+d3-d2)-q2;
fx=[P2;P3;Q2];
j11=diff(P2,d2);
j12=diff(P2,d3);
j13=diff(P2,v2);
j21=diff(P3,d2);
j22=diff(P3,d3);
j23=diff(P3,v2);
j31=diff(Q2,d2);
j32=diff(Q2,d3);
j33=diff(Q2,v2);
% Jacobian
J=[j11,j12,j13;j21,j22,j23;j31,j32,j33];

er=1;
iteration=1;
while max(er) > 1e-3 || iteration>50
    J1  = subs(J,[d2;d3;v2],X);
    fx1 = subs(fx,[d2;d3;v2],X);
    X1  = X - double(J1 \ fx1);
    er  = X - X1;
    X = X1;
    iteration = iteration + 1;
end
if iteration<50
    delta3=rad2deg(X(1,1))
    delta2 = rad2deg(X(2,1))
    v3 = X(3,1)
    disp(['results: delta3 = ',num2str(delta3),' delta2 = ',num2str(delta2), ' V3 = ',num2str(v3) ])
    disp(['The system converges in ', num2str(iteration), ' iterations.'])
else
    disp('The system did not converge!')
end