clear,clc;
numd=[1 -0.1];
dend=conv([1 -0.8],[1 -0.9]);
N_sim=60;
k=0:N_sim-1;
H=dimpulse(numd,dend,k+1);

a=0.9;
N=5;
[A1,L0]=lagd(a,N);
L(:,1)=L0;

for kk=2:N_sim;
    L(:,kk)=A1*L(:,kk-1);
end

%c1=L(1,:)*H;
%c2=L(2,:)*H;
%c3=L(3,:)*H;
%c4=L(4,:)*H;
%c5=L(5,:)*H;

for i = 1:N
    c(i,:)=L(i,:)*H;
end

H_model = zeros(1,N_sim)
%H_model=c1*L(1,:)+c2*L(2,:)+c3*L(3,:)+c4*L(4,:)+c5*L(5,:);
for i = 1:N
    H_model=H_model+c(i,:)*L(i,:)
end

figure
plot(k,H),grid
hold on
plot(k,H_model,'LineWidth',2,'Color',[.8 0 0])
%set(gca,'FontSize',12,'FontName','helvetica');
title(['Kreiranje impulsnog modela pomocu Laguereovih polinoma za: N = ',num2str(N),', a = ',num2str(a)]);
legend('podaci','model')
xlabel('Trenutci odabiranja')
ylabel('Impuslni odziv')