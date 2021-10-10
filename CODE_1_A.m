
%PROBLEM 1_A MATLAB CODE

INPUT_1_A

nele=size(coord,1);
fp = fopen('1_A_OUTPUT.txt','w');

%Initialising Element stiffness matrix and load vector matrix

Ke=zeros(3,3,nele);
Fe=zeros(3,nele);
coord_e=zeros(3,2);
for ele=1:4
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
    [ke_d,fe_d]=domain(xi,eta,coord_e,k,Q);
    fe_d
    
    %convection term
    sum_l=zeros(3,3);
    sum_r=zeros(3,1);
    if ele==1
        edgeno=1;
        for j=1:length(beta)
            [k_h,f_h]=gamma_h(beta(j),coord_e,h,Tinf,edgeno);
            sum_l=sum_l+k_h*weight(j);
            sum_r=sum_r+f_h*weight(j);
        end
        ke_h=sum_l;
        fe_h=sum_r
    end
    
    if ele==2
        edgeno=2;
        for j=1:length(beta)
            [k_h,f_h]=gamma_h(beta(j),coord_e,h,Tinf,edgeno);
            sum_l=sum_l+k_h*weight(j);
            sum_r=sum_r+f_h*weight(j);
        end
        ke_h=sum_l;
        fe_h=sum_r
    end
    
    if ele==3
        ke_h=zeros(3,3);
        fe_h=zeros(3,1)
    end
    
    if ele==4
        edgeno=1;
        for j=1:length(beta)
            [k_h,f_h]=gamma_h(beta(j),coord_e,h,Tinf,edgeno);
            sum_l=sum_l+k_h*weight(j);
            sum_r=sum_r+f_h*weight(j);
        end
        ke_h=sum_l;
        fe_h=sum_r
    end
    
    %Heat flux term
    sum_r=zeros(3,1);
    
    if (ele==1 || ele==2)
        fe_q=zeros(3,1)
    end
    if ele==3
        edgeno=2;
        for j=1:length(beta)
            f_q=gamma_q(beta(j),coord_e,qn,edgeno);
            sum_r=sum_r+f_q*weight(j);
        end
        fe_q=sum_r
    end
    
    if ele==4
        edgeno=2;
        for j=1:length(beta)
            f_q=gamma_q(beta(j),coord_e,qn,edgeno);
            sum_r=sum_r+f_q*weight(j);
        end
        fe_q=sum_r
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
K=zeros(7,7);
F=zeros(7,1);
for ele=1:4
    nd1=connect(ele,1);
    nd2=connect(ele,2);
    nd3=connect(ele,3);
    vec=[nd1 nd2 nd3];
    for i=1:3
        for j=1:3
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

for i=1:7
    for j=1:7
        fprintf(fp,'%14.4e\t',K(i,j));
    end
    fprintf(fp,'\n');
end
fprintf(fp,'\n\n=======================\n');
fprintf(fp,'Global load matrix:\n');
fprintf(fp,'=======================\n\n');

for i=1:7
   
        fprintf(fp,'%14.4e\n',F(i));
    
end
    

%Imposition of BC
active=[2 4:7];
K_reduce=K([2 4:7],[2 4:7])
F_reduce=zeros(5,1);

for i=1:length(active)
    F_reduce(i)=F(active(i))-K(active(i),1)*T1-K(active(i),3)*T1;
end
F_reduce

%Solution
T_reduce=K_reduce\F_reduce
T=zeros(7,1);
T(1)=T1; T(3)=T1;
T([2 4:7])=T_reduce

fprintf(fp,'\n\n=======================\n');
fprintf(fp,' Temperature\n ');
fprintf(fp,'=======================\n');
for i=1:7
    
fprintf(fp,'%.10d\n',T(i)); 
end
fclose(fp);


