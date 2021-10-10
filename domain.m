function [ke,fe]=domain(xi,eta,coord,k,Q)
% domain term integrand
x=coord(:,1);
y=coord(:,2);
x13=x(1)-x(3);
y13=y(1)-y(3);
x23=x(2)-x(3);
y23=y(2)-y(3);
N=[xi eta 1-xi-eta];
J=[x13 y13;x23 y23];
B=J\[1 0 -1;0 1 -1];
ke=k*(B'*B)*det(J)/2;
fe=N'*Q*det(J)/2;




