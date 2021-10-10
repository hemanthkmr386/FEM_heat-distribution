function [ke,fe]=domain(xi,eta,coord,k,Q)
% domain term integrand
%----------------------
%N vector
N=[(1-xi)*(1-eta)/4 (1+xi)*(1-eta)/4 (1+xi)*(1+eta)/4 (1-xi)*(1+eta)/4];
%dN/dxi matrix
dNdxi=0.25*[-(1-eta) 1-eta 1+eta -(1+eta);
    -(1-xi) -(1+xi) 1+xi 1-xi];
%Jacobian
J=dNdxi*coord;
%B matrix
B=J\dNdxi;
%ke and fe
ke=k*(B'*B)*det(J);
fe=N'*Q*det(J);