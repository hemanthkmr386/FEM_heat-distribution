function [ke,fe]=gamma_h(beta,coord,h,Tinf,edgeno)
x=coord(:,1);
y=coord(:,2);
xi=(beta+1)/2;
if edgeno==1
    N=[1-xi xi 0];
    l=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
    elseif edgeno==2
        N=[0 1-xi xi];
        l=sqrt((x(3)-x(2))^2+(y(3)-y(2))^2);
        elseif edgeno==3
            N=[xi 0 1-xi];
            l=sqrt((x(1)-x(3))^2+(y(1)-y(3))^2);
end
ke=h*(N'*N)*l*0.5;
fe=h*(N')*Tinf*l*0.5;
    
    
