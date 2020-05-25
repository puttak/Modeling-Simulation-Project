function F=Equations(x,Az,Bz,Ar,Br,u0,C0,Pt,epsilon,Density_bed,Deffz,Deffr,r_nodes,kg,hg,as,P,Components,RxnKinetic,R)
C_C2H6=reshape(x(1:Nz*Nr),Nz,Nr);
C_C2H4=reshape(x(Nz*Nr+1:2*Nz*Nr),Nz,Nr);
C_O2=reshape(x(2*Nz*Nr+1:3*Nz*Nr),Nz,Nr);
C_CO2=reshape(x(3*Nz*Nr+1:4*Nz*Nr),Nz,Nr);
C_CO=reshape(x(4*Nz*Nr+1:5*Nz*Nr),Nz,Nr);
C_H2O=reshape(x(5*Nz*Nr+1:6*Nz*Nr),Nz,Nr);
C_N2=reshape(x(6*Nz*Nr+1:7*Nz*Nr),Nz,Nr);
Density_fluid=reshape(x(7*Nz*Nr+1:8*Nz*Nr),Nz,Nr);
Cs_C2H6=reshape(x(8*Nz*Nr+1:9*Nz*Nr),Nz,Nr);
Cs_C2H4=reshape(x(9*Nz*Nr+1:10*Nz*Nr),Nz,Nr);
Cs_O2=reshape(x(10*Nz*Nr+1:11*Nz*Nr),Nz,Nr);
Cs_CO2=reshape(x(11*Nz*Nr+1:12*Nz*Nr),Nz,Nr);
Cs_CO=reshape(x(12*Nz*Nr+1:13*Nz*Nr),Nz,Nr);
Cs_H2O=reshape(x(13*Nz*Nr+1:14*Nz*Nr),Nz,Nr);
Cpf=reshape(x(14*Nz*Nr+1:15*Nz*Nr),Nz,Nr);
T=reshape(x(15*Nz*Nr+1:16*Nz*Nr),Nz,Nr);
Ts=reshape(x(16*Nz*Nr+1:17*Nz*Nr),Nz,Nr);

E_C_C2H6=zeros(Nz,Nr);
E_C_C2H4=zeros(Nz,Nr);
E_C_O2=zeros(Nz,Nr);
E_C_CO2=zeros(Nz,Nr);
E_C_CO=zeros(Nz,Nr);
E_C_H2O=zeros(Nz,Nr);
E_C_N2=zeros(Nz,Nr);
E_Rof=zeros(Nz,Nr);
E_Cs_C2H6=zeros(Nz,Nr);
E_Cs_C2H4=zeros(Nz,Nr);
E_Cs_O2=zeros(Nz,Nr);
E_Cs_CO2=zeros(Nz,Nr);
E_Cs_CO=zeros(Nz,Nr);
E_Cs_H2O=zeros(Nz,Nr);
E_Cpf=zeros(Nz,Nr);
E_T=zeros(Nz,Nr);
E_Ts=zeros(Nz,Nr);

for i=1:Nz
    for k=1:Nr 
        E_C_C2H6(i,k) = -u0*Az(i+1,2:end-1)*C_C2H6(:,k) - u0*Az(i+1,1)*(((u0*C0(1))+(epsilon*Deffz*Az(1,2:end-1)*C_C2H6(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C2H6(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - u0*Az(i+1,end)*(((Az(end,2:end-1)*C_C2H6(:,k))+(Az(end,1)*(((u0*C0(1))+(epsilon*Deffz*Az(1,2:end-1)*C_C2H6(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C2H6(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + epsilon*Deffz*Bz(i+1,2:end-1)*C_C2H6(:,k) + epsilon*Deffz*Bz(i+1,1)*(((u0*C0(1))+(epsilon*Deffz*Az(1,2:end-1)*C_C2H6(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C2H6(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - epsilon*Deffz*Bz(i+1,end)*(((Az(end,2:end-1)*C_C2H6(:,k))+(Az(end,1)*(((u0*C0(1))+(epsilon*Deffz*Az(1,2:end-1)*C_C2H6(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C2H6(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,2:end-1)*C_C2H6(i,:) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,1)*(((Ar(1,2:end-1)*C_C2H6(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C2H6(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,end)*(((Ar(end,2:end-1)*C_C2H6(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C2H6(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C2H6(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (epsilon*Deffr)*Br(k+1,2:end-1)*C_C2H6(i,:) + (epsilon*Deffr)*Br(k+1,1)*(((Ar(1,2:end-1)*C_C2H6(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C2H6(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + (epsilon*Deffr)*Br(k+1,end)*(((Ar(end,2:end-1)*C_C2H6(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C2H6(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C2H6(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (1-epsilon)*kg*as*(Cs_C2H6(i,k)-C_C2H6(i,k));
        
        E_C_C2H4(i,k)= -u0*Az(i+1,2:end-1)*C_C2H4(:,k) - u0*Az(i+1,1)*(((u0*C0(2))+(epsilon*Deffz*Az(1,2:end-1)*C_C2H4(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C2H4(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - u0*Az(i+1,end)*(((Az(end,2:end-1)*C_C2H4(:,k))+(Az(end,1)*(((u0*C0(2))+(epsilon*Deffz*Az(1,2:end-1)*C_C2H4(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C2H4(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + epsilon*Deffz*Bz(i+1,2:end-1)*C_C2H4(:,k) + epsilon*Deffz*Bz(i+1,1)*(((u0*C0(2))+(epsilon*Deffz*Az(1,2:end-1)*C_C2H4(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C2H4(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - epsilon*Deffz*Bz(i+1,end)*(((Az(end,2:end-1)*C_C2H4(:,k))+(Az(end,1)*(((u0*C0(2))+(epsilon*Deffz*Az(1,2:end-1)*C_C2H4(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_C2H4(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,2:end-1)*C_C2H4(i,:) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,1)*(((Ar(1,2:end-1)*C_C2H4(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C2H4(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,end)*(((Ar(end,2:end-1)*C_C2H4(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C2H4(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C2H4(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (epsilon*Deffr)*Br(k+1,2:end-1)*C_C2H4(i,:) + (epsilon*Deffr)*Br(k+1,1)*(((Ar(1,2:end-1)*C_C2H4(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C2H4(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + (epsilon*Deffr)*Br(k+1,end)*(((Ar(end,2:end-1)*C_C2H4(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_C2H4(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_C2H4(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (1-epsilon)*kg*as*(Cs_C2H4(i,k)-C_C2H4(i,k));
        
        E_C_O2(i,k)= -u0*Az(i+1,2:end-1)*C_O2(:,k) - u0*Az(i+1,1)*(((u0*C0(3))+(epsilon*Deffz*Az(1,2:end-1)*C_O2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_O2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - u0*Az(i+1,end)*(((Az(end,2:end-1)*C_O2(:,k))+(Az(end,1)*(((u0*C0(3))+(epsilon*Deffz*Az(1,2:end-1)*C_O2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_O2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + epsilon*Deffz*Bz(i+1,2:end-1)*C_O2(:,k) + epsilon*Deffz*Bz(i+1,1)*(((u0*C0(3))+(epsilon*Deffz*Az(1,2:end-1)*C_O2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_O2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - epsilon*Deffz*Bz(i+1,end)*(((Az(end,2:end-1)*C_O2(:,k))+(Az(end,1)*(((u0*C0(3))+(epsilon*Deffz*Az(1,2:end-1)*C_O2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_O2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,2:end-1)*C_O2(i,:) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,1)*(((Ar(1,2:end-1)*C_O2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_O2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,end)*(((Ar(end,2:end-1)*C_O2(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_O2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_O2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (epsilon*Deffr)*Br(k+1,2:end-1)*C_O2(i,:) + (epsilon*Deffr)*Br(k+1,1)*(((Ar(1,2:end-1)*C_O2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_O2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + (epsilon*Deffr)*Br(k+1,end)*(((Ar(end,2:end-1)*C_O2(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_O2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_O2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (1-epsilon)*kg*as*(Cs_O2(i,k)-C_O2(i,k));
        
        E_C_CO2(i,k)= -u0*Az(i+1,2:end-1)*C_CO2(:,k) - u0*Az(i+1,1)*(((u0*C0(4))+(epsilon*Deffz*Az(1,2:end-1)*C_CO2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_CO2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - u0*Az(i+1,end)*(((Az(end,2:end-1)*C_CO2(:,k))+(Az(end,1)*(((u0*C0(4))+(epsilon*Deffz*Az(1,2:end-1)*C_CO2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_CO2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + epsilon*Deffz*Bz(i+1,2:end-1)*C_CO2(:,k) + epsilon*Deffz*Bz(i+1,1)*(((u0*C0(4))+(epsilon*Deffz*Az(1,2:end-1)*C_CO2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_CO2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - epsilon*Deffz*Bz(i+1,end)*(((Az(end,2:end-1)*C_CO2(:,k))+(Az(end,1)*(((u0*C0(4))+(epsilon*Deffz*Az(1,2:end-1)*C_CO2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_CO2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,2:end-1)*C_CO2(i,:) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,1)*(((Ar(1,2:end-1)*C_CO2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_CO2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,end)*(((Ar(end,2:end-1)*C_CO2(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_CO2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_CO2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (epsilon*Deffr)*Br(k+1,2:end-1)*C_CO2(i,:) + (epsilon*Deffr)*Br(k+1,1)*(((Ar(1,2:end-1)*C_CO2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_CO2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + (epsilon*Deffr)*Br(k+1,end)*(((Ar(end,2:end-1)*C_CO2(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_CO2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_CO2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (1-epsilon)*kg*as*(Cs_CO2(i,k)-C_CO2(i,k));
        
        E_C_CO(i,k)= -u0*Az(i+1,2:end-1)*C_CO(:,k) - u0*Az(i+1,1)*(((u0*C0(5))+(epsilon*Deffz*Az(1,2:end-1)*C_CO(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_CO(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - u0*Az(i+1,end)*(((Az(end,2:end-1)*C_CO(:,k))+(Az(end,1)*(((u0*C0(5))+(epsilon*Deffz*Az(1,2:end-1)*C_CO(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_CO(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + epsilon*Deffz*Bz(i+1,2:end-1)*C_CO(:,k) + epsilon*Deffz*Bz(i+1,1)*(((u0*C0(5))+(epsilon*Deffz*Az(1,2:end-1)*C_CO(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_CO(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - epsilon*Deffz*Bz(i+1,end)*(((Az(end,2:end-1)*C_CO(:,k))+(Az(end,1)*(((u0*C0(5))+(epsilon*Deffz*Az(1,2:end-1)*C_CO(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_CO(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,2:end-1)*C_CO(i,:) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,1)*(((Ar(1,2:end-1)*C_CO(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_CO(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,end)*(((Ar(end,2:end-1)*C_CO(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_CO(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_CO(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (epsilon*Deffr)*Br(k+1,2:end-1)*C_CO(i,:) + (epsilon*Deffr)*Br(k+1,1)*(((Ar(1,2:end-1)*C_CO(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_CO(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + (epsilon*Deffr)*Br(k+1,end)*(((Ar(end,2:end-1)*C_CO(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_CO(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_CO(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (1-epsilon)*kg*as*(Cs_CO(i,k)-C_CO(i,k));
        
        E_C_H2O(i,k)= -u0*Az(i+1,2:end-1)*C_H2O(:,k) - u0*Az(i+1,1)*(((u0*C0(6))+(epsilon*Deffz*Az(1,2:end-1)*C_H2O(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_H2O(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - u0*Az(i+1,end)*(((Az(end,2:end-1)*C_H2O(:,k))+(Az(end,1)*(((u0*C0(6))+(epsilon*Deffz*Az(1,2:end-1)*C_H2O(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_H2O(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + epsilon*Deffz*Bz(i+1,2:end-1)*C_H2O(:,k) + epsilon*Deffz*Bz(i+1,1)*(((u0*C0(6))+(epsilon*Deffz*Az(1,2:end-1)*C_H2O(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_H2O(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - epsilon*Deffz*Bz(i+1,end)*(((Az(end,2:end-1)*C_H2O(:,k))+(Az(end,1)*(((u0*C0(6))+(epsilon*Deffz*Az(1,2:end-1)*C_H2O(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_H2O(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,2:end-1)*C_H2O(i,:) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,1)*(((Ar(1,2:end-1)*C_H2O(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_H2O(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,end)*(((Ar(end,2:end-1)*C_H2O(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_H2O(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_H2O(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (epsilon*Deffr)*Br(k+1,2:end-1)*C_H2O(i,:) + (epsilon*Deffr)*Br(k+1,1)*(((Ar(1,2:end-1)*C_H2O(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_H2O(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + (epsilon*Deffr)*Br(k+1,end)*(((Ar(end,2:end-1)*C_H2O(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_H2O(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_H2O(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (1-epsilon)*kg*as*(Cs_H2O(i,k)-C_H2O(i,k));
        
        E_C_N2(i,k)= -u0*Az(i+1,2:end-1)*C_N2(:,k) - u0*Az(i+1,1)*(((u0*C0(7))+(epsilon*Deffz*Az(1,2:end-1)*C_N2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_N2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - u0*Az(i+1,end)*(((Az(end,2:end-1)*C_N2(:,k))+(Az(end,1)*(((u0*C0(7))+(epsilon*Deffz*Az(1,2:end-1)*C_N2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_N2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + epsilon*Deffz*Bz(i+1,2:end-1)*C_N2(:,k) + epsilon*Deffz*Bz(i+1,1)*(((u0*C0(7))+(epsilon*Deffz*Az(1,2:end-1)*C_N2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_N2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1)))) - epsilon*Deffz*Bz(i+1,end)*(((Az(end,2:end-1)*C_N2(:,k))+(Az(end,1)*(((u0*C0(7))+(epsilon*Deffz*Az(1,2:end-1)*C_N2(:,k))-(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,2:end-1)*C_N2(:,k)))/((u0)+(((epsilon*Deffz*Az(1,end))/(Az(end,end)))*Az(end,1))-(epsilon*Deffz*Az(1,1))))))/(-Az(end,end))) ...
            + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,2:end-1)*C_N2(i,:) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,1)*(((Ar(1,2:end-1)*C_N2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_N2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + ((epsilon*Deffr)/r_nodes(k+1))*Ar(k+1,end)*(((Ar(end,2:end-1)*C_N2(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_N2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_N2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end))) ...
            + (epsilon*Deffr)*Br(k+1,2:end-1)*C_N2(i,:) + (epsilon*Deffr)*Br(k+1,1)*(((Ar(1,2:end-1)*C_N2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_N2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1)))) + (epsilon*Deffr)*Br(k+1,end)*(((Ar(end,2:end-1)*C_N2(i,:))+(Ar(end,1)*(((Ar(1,2:end-1)*C_N2(i,:))-((Ar(1,end)/Ar(end,end))*Ar(end,2:end-1)*C_N2(i,:)))/(((Ar(1,end)/Ar(end,end))*Ar(end,1))-(Ar(1,1))))))/(-Ar(end,end)));
        
        E_Rof(i,k)= Density_fluid(i,k) - (P*((C_C2H6(i,k)/(C_C2H6(i,k)+C_C2H4(i,k)+C_O2(i,k)+C_CO2(i,k)+C_CO(i,k)+C_H2O(i,k)+C_N2(i,k)))*Components(1).Mw + (C_C2H4(i,k)/(C_C2H6(i,k)+C_C2H4(i,k)+C_O2(i,k)+C_CO2(i,k)+C_CO(i,k)+C_H2O(i,k)+C_N2(i,k)))*Components(2).Mw + (C_O2(i,k)/(C_C2H6(i,k)+C_C2H4(i,k)+C_O2(i,k)+C_CO2(i,k)+C_CO(i,k)+C_H2O(i,k)+C_N2(i,k)))*Components(3).Mw + (C_CO2(i,k)/(C_C2H6(i,k)+C_C2H4(i,k)+C_O2(i,k)+C_CO2(i,k)+C_CO(i,k)+C_H2O(i,k)+C_N2(i,k)))*Components(4).Mw + (C_CO(i,k)/(C_C2H6(i,k)+C_C2H4(i,k)+C_O2(i,k)+C_CO2(i,k)+C_CO(i,k)+C_H2O(i,k)+C_N2(i,k)))*Components(5).Mw ...
            + (C_H2O(i,k)/(C_C2H6(i,k)+C_C2H4(i,k)+C_O2(i,k)+C_CO2(i,k)+C_CO(i,k)+C_H2O(i,k)+C_N2(i,k)))*Components(6).Mw + (C_N2(i,k)/(C_C2H6(i,k)+C_C2H4(i,k)+C_O2(i,k)+C_CO2(i,k)+C_CO(i,k)+C_H2O(i,k)+C_N2(i,k)))*Components(7).Mw))/(R*T(i,k)); % [g/m^3] Density of fluid
        
        C_gas   = [C_C2H6(i,k) C_C2H4(i,k) C_O2(i,k) C_CO2(i,k) C_CO(i,k) C_H2O(i,k) C_N2(i,k)];
        y_gas   = C_gas/sum(C_gas);   % mole fraction of components in order: [C2H6 C2H4 O2 CO2 CO H2O N2]
        
        C_solid = [Cs_C2H6(i,k) Cs_C2H4(i,k) Cs_O2(i,k) Cs_CO2(i,k) Cs_CO(i,k) Cs_H2O(i,k) Cs_N2(i,k)];
        
        E_Cs_C2H6(i,k)= (1-epsilon)*kg*as*(C_C2H6(i,k) - Cs_C2H6(i,k)) + Density_bed*(ODHReactions(C_solid,Ts(i,k),R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],1));
        
        E_Cs_C2H4(i,k)= (1-epsilon)*kg*as*(C_C2H6(i,k) - Cs_C2H6(i,k)) + Density_bed*(ODHReactions(C_solid,Ts(i,k),R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],2));
        
        E_Cs_O2(i,k)=   (1-epsilon)*kg*as*(C_C2H6(i,k) - Cs_C2H6(i,k)) + Density_bed*(ODHReactions(C_solid,Ts(i,k),R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],3));
        
        E_Cs_CO2(i,k)=  (1-epsilon)*kg*as*(C_C2H6(i,k) - Cs_C2H6(i,k)) + Density_bed*(ODHReactions(C_solid,Ts(i,k),R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],4));
        
        E_Cs_CO(i,k)=   (1-epsilon)*kg*as*(C_C2H6(i,k) - Cs_C2H6(i,k)) + Density_bed*(ODHReactions(C_solid,Ts(i,k),R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],5));
        
        E_Cs_H2O(i,k)=  (1-epsilon)*kg*as*(C_C2H6(i,k) - Cs_C2H6(i,k)) + Density_bed*(ODHReactions(C_solid,Ts(i,k),R,Pt,Flowin,RxnKinetic,[Components.deltaS0],[Components.deltaH0],6));
        
        E_C_N2(i,k)=    (1-epsilon)*kg*as*(C_C2H6(i,k) - Cs_C2H6(i,k)); % Term of reaction must be ommited for N2 component, consider it's an inert component
        
        E_Cpf= Cpf(i,k) - (((Components(1).cp_R(1) + Components(1).cp_R(2)*T(i,k) + Components(1).cp_R(3)*T(i,k)^2 + Components(1).cp_R(4)*T(i,k)^(-2))*R)*(y_gas(1))) ...
            - (((Components(2).cp_R(1) + Components(2).cp_R(2)*T(i,k) + Components(2).cp_R(3)*T(i,k)^2 + Components(2).cp_R(4)*T(i,k)^(-2))*R)*(y_gas(2))) ...
            - (((Components(3).cp_R(1) + Components(3).cp_R(2)*T(i,k) + Components(3).cp_R(3)*T(i,k)^2 + Components(3).cp_R(4)*T(i,k)^(-2))*R)*(y_gas(3))) ...
            - (((Components(4).cp_R(1) + Components(4).cp_R(2)*T(i,k) + Components(4).cp_R(3)*T(i,k)^2 + Components(4).cp_R(4)*T(i,k)^(-2))*R)*(y_gas(4))) ...
            - (((Components(5).cp_R(1) + Components(5).cp_R(2)*T(i,k) + Components(5).cp_R(3)*T(i,k)^2 + Components(5).cp_R(4)*T(i,k)^(-2))*R)*(y_gas(5))) ...
            - (((Components(6).cp_R(1) + Components(6).cp_R(2)*T(i,k) + Components(6).cp_R(3)*T(i,k)^2 + Components(6).cp_R(4)*T(i,k)^(-2))*R)*(y_gas(6))) ...
            - (((Components(7).cp_R(1) + Components(7).cp_R(2)*T(i,k) + Components(7).cp_R(3)*T(i,k)^2 + Components(7).cp_R(4)*T(i,k)^(-2))*R)*(y_gas(7)));
        
        E_T(i,k)= ;
        
        E_Ts(i,k)= (1-epsilon)*hg*as*(T(i,k) - Ts(i,k)) + Density_bed;
    end
end

F=[reshape(E_C_C2H6,Nz*Nr,1)  ;  reshape(E_C_C2H4,Nz*Nr,1)  ; ...
   reshape(E_C_O2,Nz*Nr,1)    ;  reshape(E_C_CO2,Nz*Nr,1)   ; ...
   reshape(E_C_CO,Nz*Nr,1)    ;  reshape(E_C_H2O,Nz*Nr,1)   ; ...
   reshape(E_C_N2,Nz*Nr,1)    ;  reshape(E_Rof,Nz*Nr,1)     ; ...
   reshape(E_Cs_C2H6,Nz*Nr,1) ;  reshape(E_Cs_C2H4,Nz*Nr,1) ; ...
   reshape(E_Cs_O2,Nz*Nr,1)   ;  reshape(E_Cs_CO2,Nz*Nr,1)  ; ...
   reshape(E_Cs_CO,Nz*Nr,1)   ;  reshape(E_Cs_H2O,Nz*Nr,1)  ; ...
   reshape(E_Cpf,Nz*Nr,1)     ;reshape(E_T,Nz*Nr,1)         ; ... 
   reshape(E_Ts,Nz*Nr,1)     ];
  















