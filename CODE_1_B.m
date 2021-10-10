
%MATLAB CODE FOR PROBLEM 1_B

INPUT_1_B

nele=size(connect,1);
fp = fopen('1_B_OUTPUT.txt','w');

%Initialising Element stiffness matrix and load vector matrix

Ke=zeros(6,6,nele);
Fe=zeros(6,nele);
coord_e=zeros(6,2);

for ele=1:nele
    
    %condition for heat generation
    if (ele==1 || ele==2)
        Q=Q1;
    else
        Q=Q2;
    end
    
    % nodes of each element
    node=connect(ele,:)
    % coordinates of each element
    for j=1:length(node)
        coord_e(j,:)=coord(node(j),:);
    end
    coord_e
    % domain term
    sum_l=zeros(6,6);
    sum_r=zeros(6,1);
    for i=1:length(xi)
        [k_d,f_d]=domain(xi(i),eta(i),coord_e,k,Q);
        sum_l=sum_l+k_d*w(i);
        sum_r=sum_r+f_d*w(i);
    end
    ke_d=sum_l;
    fe_d=sum_r;
    
    %convection term
    sum_l=zeros(6,6); %variable to calculate Gaussian quadrature
    sum_r=zeros(6,1);
    if ele==1
        edgeno=1;
        for j=1:length(beta)
            [k_h,f_h]=gamma_h(beta(j),coord_e,h,Tinf,edgeno);
            sum_l=sum_l+k_h*weight(j);
            sum_r=sum_r+f_h*weight(j);
        end
        ke_h=sum_l;
        fe_h=sum_r;
    end
    
    if ele==2
        edgeno=2;
        for j=1:length(beta)
            [k_h,f_h]=gamma_h(beta(j),coord_e,h,Tinf,edgeno);
            sum_l=sum_l+k_h*weight(j);
            sum_r=sum_r+f_h*weight(j);
        end
        ke_h=sum_l;
        fe_h=sum_r;
    end
    
    if ele==3
        ke_h=zeros(6,6);
        fe_h=zeros(6,1);
    end
    
    if ele==4
        edgeno=1;
        for j=1:length(beta)
            [k_h,f_h]=gamma_h(beta(j),coord_e,h,Tinf,edgeno);
            sum_l=sum_l+k_h*weight(j);
            sum_r=sum_r+f_h*weight(j);
        end
        ke_h=sum_l;
        fe_h=sum_r;
    end
    
    %Heat flux term
    sum_r=zeros(6,1);
    
    if (ele==1 || ele==2)
        fe_q=zeros(6,1);
    end
    if ele==3
        edgeno=2;
        for j=1:length(beta)
            f_q=gamma_q(beta(j),coord_e,qn,edgeno);
            sum_r=sum_r+f_q*weight(j);
        end
        fe_q=sum_r;
    end
    
    if ele==4
        edgeno=2;
        for j=1:length(beta)
            f_q=gamma_q(beta(j),coord_e,qn,edgeno);
            sum_r=sum_r+f_q*weight(j);
        end
        fe_q=sum_r;
    end
    %Element stiffness matrix
    Ke(:,:,ele)=Ke(:,:,ele)+ke_d+ke_h;
    Fe(:,ele)=Fe(:,ele)+fe_d+fe_h+fe_q;
end

k1=Ke(:,:,1)
k2=Ke(:,:,2)
k3=Ke(:,:,3)
k4=Ke(:,:,4)

f1=Fe(:,1)
f2=Fe(:,2)
f3=Fe(:,3)
f4=Fe(:,4)

% Assembly
K=zeros(16,16);
F=zeros(16,1);
for ele=1:nele
    nd1=connect(ele,1);
    nd2=connect(ele,2);
    nd3=connect(ele,3);
    nd4=connect(ele,4);
    nd5=connect(ele,5);
    nd6=connect(ele,6);
    vec=[nd1 nd2 nd3 nd4 nd5 nd6];
    for i=1:6
        for j=1:6
            K(vec(i),vec(j))=K(vec(i),vec(j))+Ke(i,j,ele);
        end
        F(vec(i))=F(vec(i))+Fe(i,ele);
    end
end
K
F

fprintf(fp,'\n\n=======================\n');
fprintf(fp,'Global stiffness matrix:\n');
fprintf(fp,'=======================\n\n');

for i=1:16
    for j=1:16
        fprintf(fp,'%14.4e\t',K(i,j));
    end
    fprintf(fp,'\n');
end
fprintf(fp,'\n\n=======================\n');
fprintf(fp,'Global load matrix:\n');
fprintf(fp,'=======================\n\n');

for i=1:16
   
        fprintf(fp,'%14.4e\n',F(i));
    
end

%Imposition of BC
active=[2 4:9 11:16]; %active DOF
K_reduce=K([2 4:9 11:16],[2 4:9 11:16])
F_reduce=zeros(13,1);

for i=1:length(active)
    F_reduce(i)=F(active(i))-K(active(i),1)*T1-K(active(i),3)*T1-K(active(i),10)*T1;
end
F_reduce

%Solution
T_reduce=K_reduce\F_reduce
T=zeros(16,1);
T(1)=T1; T(3)=T1; T(10)=T1;
T([2 4:9 11:16])=T_reduce;

fprintf(fp,'\n\n=======================\n');
fprintf(fp,' Temperature\n ');
fprintf(fp,'=======================\n');
for i=1:16
    
fprintf(fp,'%.10d\n',T(i)); 
end
fclose(fp);



