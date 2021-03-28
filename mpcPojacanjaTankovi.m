clear, clc;
 
Ac = [-1 1;
     1  -2];
 
Bc = [1 0;
     0 1];
 
Cc = [1 0;
     0 1];

Dc = [0 0;
     0 0];
 
 
sis = ss(Ac,Bc,Cc,Dc)

sisd = c2d(sis,0.1)
[Ap,Bp,Cp,Dp] = ssdata(sisd)
 
Np=5;            %Horizont predikcije stanja (izlaza)
Nc=4;            %Horizont predikcije upravljanja
r_omega=0.1;     %Tezinski koeficient upravljanja

[m1,n1] = size(Cp);
[n1,n_in] = size(Bp);

%Kreiranej izmenjenih matrica modela u odnosu na odstupanja
A_e = eye(n1+m1,n1+m1);
A_e(1:n1,1:n1) = Ap;
A_e(n1+1:n1+m1,1:n1) = Cp*Ap;

B_e = zeros(n1+m1,n_in);
B_e(1:n1,:) = Bp;
B_e(n1+1:n1+m1,:) = Cp*Bp;

C_e = zeros(m1,n1+m1);
C_e(:,n1+1:n1+m1) = eye(m1,m1);

[m1,n1] = size(C_e);
[n1,n_in] = size(B_e);

n = n1+m1;

% Kreiranje matrice F = [CA;CA^2;...;CA^Np]
F = zeros(m1*Np,n1);
stepen = 1
for kk = 1:m1:Np*m1   
    F(kk:kk+m1-1,:) = C_e*A_e^stepen;
    stepen = stepen +1;
end

% Kreiranje matrice Phi koja je Toeplitz-ova matrica
h = zeros(m1*Nc,n1);
%h(1:m1,:) = C_e
stepen = 0;
for kk = 1:m1:Np*m1   
    h(kk:kk+m1-1,:) = C_e*A_e^stepen;
    stepen = stepen +1;
end

v = h*B_e;
Phi = zeros(Np*m1,Nc*m1); % dimezija Phi
Phi(:,1:m1) = v(:,1:m1); %Prva kolona Phi


for i = m1+1:n_in:Nc*n_in
    Phi(:,i:i+m1-1)=[zeros(i-1,m1);v(1:Np*m1-i+1,:)]; %Toeplitz-ova matrica
end



% Kreiranje ostalih matrica za predikciju
%BarRs - Rs nadvuceno
BarRs = zeros(Np*m1,n_in);
for i = 1:m1:Np*m1
    BarRs(i:i+m1-1,:) = eye(n_in);
end
Phi_Phi = Phi'*Phi;
Phi_F = Phi'*F;  
Phi_R = Phi'*BarRs;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Provera
%Fmoje = [C_e*A_e; C_e*A_e^2; C_e*A_e^3; C_e*A_e^4; C_e*A_e^5;]

%o = [0 0;0 0];

%Fi = [   C_e*B_e          o                  o                  o               o;
%      C_e*A_e^1*B_e     C_e*B_e             o                  o                o;
%      C_e*A_e^2*B_e     C_e*A_e^1*B_e      C_e*B_e             o                o;
%      C_e*A_e^3*B_e     C_e*A_e^2*B_e      C_e*A_e^1*B_e     C_e*B_e            o;
%      C_e*A_e^4*B_e     C_e*A_e^3*B_e      C_e*A_e^2*B_e     C_e*A_e^1*B_e     C_e*B_e];














