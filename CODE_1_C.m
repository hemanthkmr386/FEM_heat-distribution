
%MATLAB CODE FOR PROBLEM 1_C

INPUT_1_C

nele=size(connect,1);

nnode=size(connect,2);

fp = fopen('1_C_OUTPUT.txt','w');


%Initialising with zeros

%ELEMENT 1

coord_1=zeros(nnode,2);
node1=connect(1,:);
for i=1:nnode
    coord_1(i,:)=coord(node1(i),:);
end

%Initialising with zeros

%DOMAIN TERM
sum_l=zeros(nnode,nnode);
sum_r=zeros(nnode,1);
for i=1:length(xi)
    for j=1:length(eta)
        [k1_d,f1_d]=domain(xi(i),eta(j),coord_1,k,Q1);
        sum_l=sum_l+k1_d*xiw(i)*etaw(j);
        sum_r=sum_r+f1_d*xiw(i)*etaw(j);
    end
end
ke1_d=sum_l;
fe1_d=sum_r;

%CONVECTION TERM
%Initialising with zeros

sum_l=zeros(nnode,nnode);
sum_r=zeros(nnode,1);
for i=1:length(xi)
    [k1_h,f1_h]=gamma_h(xi(i),coord_1,h,Tinf,1);
    sum_l=sum_l+k1_h*xiw(i);
    sum_r=sum_r+f1_h*xiw(i);
end
k1_h1=sum_l;
f1_h1=sum_r;

sum_l=zeros(nnode,nnode);
sum_r=zeros(nnode,1);
for i=1:length(xi)
    [k1_h,f1_h]=gamma_h(xi(i),coord_1,h,Tinf,3);
    sum_l=sum_l+k1_h*xiw(i);
    sum_r=sum_r+f1_h*xiw(i);
end
k1_h2=sum_l;
f1_h2=sum_r;

ke1_h=k1_h1+k1_h2;
fe1_h=f1_h1+f1_h2;

ke1=ke1_d+ke1_h
fe1=fe1_d+fe1_h

%ELEMENT 2
%-------------------------
coord_2=zeros(nnode,2);
node2=connect(2,:);
for i=1:nnode
    coord_2(i,:)=coord(node2(i),:);
end

%DOMAIN TERM
sum_l=zeros(nnode,nnode);
sum_r=zeros(nnode,1);
for i=1:length(xi)
    for j=1:length(eta)
        [k2_d,f2_d]=domain(xi(i),eta(j),coord_2,k,Q2);
        sum_l=sum_l+k2_d*xiw(i)*etaw(j);
        sum_r=sum_r+f2_d*xiw(i)*etaw(j);
    end
end
ke2_d=sum_l;
fe2_d=sum_r;

%CONVECTION TERM
sum_l=zeros(nnode,nnode);
sum_r=zeros(nnode,1);
for i=1:length(xi)
    [k2_h,f2_h]=gamma_h(xi(i),coord_2,h,Tinf,1);
    sum_l=sum_l+k2_h*xiw(i);
    sum_r=sum_r+f2_h*xiw(i);
end
ke2_h=sum_l;
fe2_h=sum_r;

%HEAT FLUX TERM
sum_r=zeros(nnode,1);
for i=1:length(xi)
    f2_q=gamma_q(xi(i),coord_2,qn,2);
    sum_r=sum_r+f2_q*xiw(i);
end
f2_q1=sum_r;

sum_r=zeros(nnode,1);
for i=1:length(xi)
    f2_q=gamma_q(xi(i),coord_2,qn,3);
    sum_r=sum_r+f2_q*xiw(i);
end
f2_q2=sum_r;

fe2_q=f2_q1+f2_q2;

ke2=ke2_d+ke2_h
fe2=fe2_d+fe2_h+fe2_q


%------------------------------------
% Assembly
%------------------------------------
K=zeros(7,7);
F=zeros(7,1);
K(1:4,1:4)=K(1:4,1:4)+ke1;
K([2,7,6,5],[2,7,6,5])=K([2,7,6,5],[2,7,6,5])+ke2

F(1:4)=F(1:4)+fe1;
F([2,7,6,5])=F([2,7,6,5])+fe2
%-------------------------------------

%-------------------------------------
%Imposition of BC
%-------------------------------------
active=[2 3 5 6 7]; %active DOF
K_reduce=K([2 3 5 6 7],[2 3 5 6 7])
F_reduce=zeros(5,1);

for i=1:length(active)
    F_reduce(i)=F(active(i))-K(active(i),1)*T1-K(active(i),4)*T1;
end
F_reduce
%------------------------------------

%Solution
T_reduce=K_reduce\F_reduce
T=zeros(7,1);
T(1)=T1; T(4)=T1;
T([2 3 5 6 7])=T_reduce

%% Writing the results in an output file
fprintf(fp,"============================================================\n");
fprintf(fp,"Element stiffness matrix for first element:\n");
fprintf(fp,"============================================================\n");
for i=1:4
    for j=1:4
        fprintf(fp,"%10.4e \t",ke1(i,j));
    end
    fprintf(fp,"\n\n");
end

fprintf(fp,"===========================================================\n");
fprintf(fp,"Element stiffness matrix for second element:\n");
fprintf(fp,"===========================================================\n");
for i=1:4
    for j=1:4
        fprintf(fp,"%10.4e \t",ke2(i,j));
    end
    fprintf(fp,"\n\n");
end

fprintf(fp,"=====================================================\n");
fprintf(fp,"Element load vector for first element:\n");
fprintf(fp,"=====================================================\n");
for i=1:4
    fprintf(fp,"%10.4e\n\n",fe1(i));
end

fprintf(fp,"=====================================================\n");
fprintf(fp,"Element load vector for second element:\n");
fprintf(fp,"=====================================================\n");
for i=1:4
    fprintf(fp,"%10.4e\n\n",fe2(i));
end

fprintf(fp,"========================================================\n");
fprintf(fp,"Global stiffness matrix:\n");
fprintf(fp,"========================================================\n");
for i=1:7
    for j=1:7
        fprintf(fp,"%10.4e \t",K(i,j));
    end
    fprintf(fp,"\n\n");
end

fprintf(fp,"=====================================================\n");
fprintf(fp,"Global load vector:\n");
fprintf(fp,"=====================================================\n");
for i=1:7
    fprintf(fp,"%10.4e\n\n",F(i));
end

fprintf(fp,"=====================================================\n");
fprintf(fp,"Temperature at the nodes:\n");
fprintf(fp,"=====================================================\n");
for i=1:7
    fprintf(fp,"%10.4e\n\n",T(i));
end
fclose(fp);
