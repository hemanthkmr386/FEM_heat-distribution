function [ke,fe]=domain(xi,eta,coord,k,Q)
% domain term integrand
alpha=1-xi-eta;
%N vector
N=[xi*(2*xi-1) eta*(2*eta-1) alpha*(2*alpha-1) 4*xi*eta 4*eta*alpha 4*alpha*xi];
%dN/dxi matrix
dNdxi=[4*xi-1 0 1-4*alpha 4*eta -4*eta 4-8*xi-4*eta;
    0 4*eta-1 1-4*alpha 4*xi 4-4*xi-8*eta -4*xi];
%Jacobian
J=dNdxi*coord;
%B matrix
B=J\dNdxi;
%ke and fe
ke=k*(B'*B)*det(J)/2;
fe=N'*Q*det(J)/2;