clear, clc;

dt = 0.01;
numc = [10];
denc = [1 0.1 3];

Gs=tf(numc,denc)

[Ac,Bc,Cc,Dc] = tf2ss(numc,denc)

[Ap,Bp,Cp,Dp] = c2dm(Ac,Bc,Cc,Dc,dt,'zoh')

Np=30; % Horizont predvidjanja
Nc=10; % Horizont upravljanja
r_omega = 0.01;
N_sim=60; % Broj trenutaka odabiranja

% Ogranicenja
u_max = 4;
u_min = -4;
deltau_max = -1.5;
deltau_min = 1.5;


[Rbar, Phi, F, Phi_Phi, Phi_F,Phi_R, A_e, B_e, C_e] = mpcPojacanjaMIMO(Ap,Bp,Cp,Nc,Np,r_omega);
[n,n_in]=size(B_e);

%Pocetne vrednosti
xm=[0;0];
Xf=zeros(n,1);

%Referentne vrednosti
r=ones(N_sim,1);
u=0; % u(k-1) =0
y=0;

% Matrice za  optimizaciju sa ogranicenjima
A_cons = zeros(4*Nc,Nc);
for i=1:Nc
    for j=1:Nc     
        if i==j
            A_cons(i,j)=1;   
            A_cons(i+Nc,j)=-1;
        end;
        if i>=j 
            A_cons(i+2*Nc,j)=1;
            A_cons(i+3*Nc,j)=-1;
        end;        
    end;
end;

b = zeros(4*Nc,1);
for i = 1:Nc
    b(i,1) = deltau_min;
    b(i+Nc,1) = -deltau_max;
    b(i+2*Nc,1) = u_max-u;
    b(i+3*Nc,1) = -u_min+u;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%Simulacija promene ulaza

for j=1:N_sim;
DeltaU1 = inv(Phi_Phi+r_omega*eye(Nc,Nc)) * (Phi_R*r(j)-Phi_F*Xf)
H = Phi_Phi+r_omega*eye(Nc,Nc);
f = -Phi_R*r(j)+Phi_F*Xf;

% Matrice za  optimizaciju sa ogranicenjima
for i = 1:Nc
    
    b(i+2*Nc,1) = u_max-u;
    b(i+3*Nc,1) = -u_min+u;
end

%A_cons = [1 0 0;
%          0 1 0;
%          0 0 1;
%         -1 0 0;
%          0 -1 0;
%          0 0 -1;
%          1 0 0;
%          1 1 0;
%          1 1 1;
%         -1 0 0;
%         -1 -1 0;
%         -1 -1 -1];

% 1.5 <= delatau(k) <= 3.0 i -3 <= u(k) <= 6 
%b = [deltau_min;
%     deltau_min;
%     deltau_min;
%     -deltau_max;
%    -deltau_max;
%    -deltau_max;
%     u_max-u;
%     u_max-u;
%     u_max-u;
%     -u_min+u;
%     -u_min+u;
%     -u_min+u];
 
DeltaU = QPhild(H,f,A_cons,b)
deltau=DeltaU(1,1);
u=u+deltau;
u1(j)=u;
y1(j)=y;
xm_old=xm;
xm=Ap*xm+Bp*u;
y=Cp*xm;
Xf=[xm-xm_old;y];
deltau_crtanje(j)=deltau;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%Crtanje grafika

k=0:(N_sim-1);

figure(1)
subplot(311)
plot(k,y1),grid
xlabel('Trenutci odabiranja')
legend('Izlaz - y')
subplot(312)
stairs(k,u1),grid
xlabel('Trenutci odabiranja')
legend('Upravljanje - u')
subplot(313)
plot(k,deltau_crtanje'),grid
xlabel('Trenutci odabiranja')
legend('Prirataj upravljanja - \Deltau')
