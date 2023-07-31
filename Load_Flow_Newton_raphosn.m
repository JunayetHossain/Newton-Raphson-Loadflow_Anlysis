clc
tic
                          %% Read Excel File Datas %%
C = xlsread('Test.xlsx'); %Reading Line Data
Load = xlsread('Test_Load.xlsx');% Reading Load Flow Data
Load(isnan(Load))=0; %Converting Nan Elements to zero
                          %% Forming Y bus Matrix %%
[row,~] = size(C);
Z = zeros(max(C(:,2)));
for j = 1:row %%Creating Impedance Matrix
    Z(C(j,1),C(j,2)) = C(j,3)+ (1i*C(j,4));
    Z(C(j,2),C(j,1)) = Z(C(j,1),C(j,2));
end
fprintf('Impedance Matrix Z =\n')
disp(Z)
y = 1./Z; %%Creating Admittance Matrix
y(y==inf)=0; %Converting infinity to zero
Y = zeros(size(y));
for m = 1:length(y) %Creating Y-bus Matrix
    sum = 0;
    for light = 1:length(y)
            sum = sum+y(m,light);     
            Y(m,light) = -y(m,light);
    end
    Y(m,m) = sum;
end
fprintf('Ybus MAatrix Y=\n')
disp(Y) %%Displaying Y-bus Matrix
              %% Seperating and defining All necessary Values %%
V = Load(:,2); V(V==0)=1; No_of_bus=length(Y);del=zeros(No_of_bus,1);
Pg = Load(:,3);Qg = Load(:,4);
PL = Load(:,5);QL = Load(:,6);
type=zeros(No_of_bus,1);
for i = 1:No_of_bus
    if i == 1
        type(i)=1;
    elseif Pg(i)>0
            type(i)=2;
    else
            type(i)=3;
    end
end
PV_bus=find(type==1|type==2);PQ_bus=find(type==3);
No_of_PVbus=length(PV_bus);No_of_PQbus=length(PQ_bus);
active_power=Pg-PL;reactive_power=Qg-QL;
iteration = 1;error = 1;
              %........Newton Raphson Load Flow Iteration........%
%=========================================================================%
while error>1e-5%Iterate until error goes below Tolerance Value(0.00001)
                          %% Forming |PQ| Matrix
    P=zeros(No_of_bus,1);
    Q=zeros(No_of_bus,1);
    for i=1:No_of_bus
        for j=1:No_of_bus
            theta = angle(Y(i,j))-del(i)+del(j);
            P(i)=P(i)+V(i)*V(j)*abs(Y(i,j))*cos(theta);
            Q(i)=Q(i)-V(i)*V(j)*abs(Y(i,j))*sin(theta);
        end
    end
    dPa=active_power-P;
    dQa=reactive_power-Q;
    dP=dPa(2:No_of_bus);
    dQ=zeros(No_of_PQbus,1);
    k=1;
    for i = 1:No_of_bus
        if type(i)==3
            dQ(k,1)=dQa(i);
            k=k+1;
        end
    end
    PQ=[dP;dQ];
                      %% Forming JACOBIAN MATRIX[H L;M N] %%
   %======================================================================%                   
   %% Formation of H(J1) 
   H=[No_of_bus-1,No_of_bus-1];
    for i = 1:No_of_bus-1
        m=i+1;
        for j = 1:No_of_bus-1
            n = j+1;
            if m == n
                for n = 1:No_of_bus
                    if m ~=n
                        theta = sin(angle(Y(m,n))-del(m)+del(n));
                        H(i,j) = H(i,j) + V(m)*V(n)*abs(Y(m,n))*theta;
                    end
                end
            else
                theta = sin(angle(Y(m,n))-del(m)+del(n));
                H(i,j) = -V(m)*V(n)*abs(Y(m,n))*theta;
            end
        end
    end
    %% Formation of L(J2)
    L=zeros(No_of_bus-1,No_of_PQbus);
    for i = 1:No_of_bus-1
        m = i+1;
        for j = 1:No_of_PQbus
            n = PQ_bus(j);
            if m~=n
                theta =cos(angle(Y(m,n))-del(m)+del(n));
                L(i,j) = V(m)*abs(Y(m,n))*theta;
            else
                for n = 1:No_of_bus
                    if m~=n
                        theta = cos(angle(Y(m,n))-del(m)+del(n));
                        L(i,j) = L(i,j) + V(n)*abs(Y(m,n))*theta;
                    end
                end
                L(i,j)=L(i,j)+2*V(m)*abs(Y(m,m))*cos(angle(Y(m,m)));
            end
        end
    end
    %% Formation of M(J3)
    M = zeros(No_of_PQbus,No_of_bus-1);
    for i = 1:No_of_PQbus
        m = PQ_bus(i);
        for j = 1:No_of_bus-1
            n = j+1;
            if m==n
                for n = 1:No_of_bus
                    if m~=n
                        theta = cos(angle(Y(m,n))-del(m)+del(n));
                        M(i,j)= M(i,j)+ V(m)*V(n)*abs(Y(m,n))*theta;
                    end
                end
            else
                theta = cos(angle(Y(m,n))-del(m)+del(n));
                M(i,j)=-V(m)*V(n)*abs(Y(m,n))*theta;
            end 
        end
    end
   %% Formation of N(J4)
   N = zeros(No_of_PQbus,No_of_PQbus);
   for i = 1:No_of_PQbus
       m = PQ_bus(i);
       for j = 1:No_of_PQbus
           n = PQ_bus(j);
           if m == n
               for n = 1:No_of_bus
                   if m~=n
                       theta = sin(angle(Y(m,n))-del(m)+del(n));
                       N(i,j)=N(i,j)+ V(n)*abs(Y(m,n))*theta;
                   end
               end
               N(i,j)=-2*V(m)*abs(Y(m,m))*sin(angle(Y(m,m)))-N(i,j);
           else
               theta = sin(angle(Y(m,n))-del(m)+del(n));
               N(i,j)=-V(m)*abs(Y(m,n))*theta;
           end
       end
   end
   J = [H L;M N]; %Formation of Jacobian Matrix
   X = J\PQ;
   dTh = X(1:No_of_bus-1);
   dV = X(No_of_bus:end);
   del(2:No_of_bus) = del(2:No_of_bus) + dTh;
   k = 1;
   for i = 1:No_of_bus
       if type(i)==3
           V(i)=V(i)+dV(k);
           k=k+1;
       end
   end
   iteration = iteration + 1;
   error = max(abs(PQ));
   
end
del = rad2deg(del);
 %% Load Flow Solution 
disp('----------------------------------------');
disp('  Newton Raphson Loadflow Solution    ');
disp('----------------------------------------');
disp(' |Bus |   |Voltage|    |Angle |');
disp(' | No.|   |pu     |    |Degree|');
disp('----------------------------------------');
for m=1:No_of_bus 
    fprintf(' %3g   ' ,m);
    fprintf(' %8.6f    ' ,V(m));
    fprintf(' %8.6f  ' ,del(m));
    fprintf('\n');
end
disp('----------------------------------------');
fprintf( 'Number Of Ieration %3g \n',iteration)
fprintf('Error %.6f\n',error);
toc      