function fe=gamma_q(beta,coord,qn,edgeno)
x=coord(:,1);
y=coord(:,2);
xi=(beta+1)/2;
if edgeno==1
    N=[(1-xi)*(1-2*xi) xi*(2*xi-1) 0 4*xi*(1-xi) 0 0];
    l=sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
    elseif edgeno==2
        N=[0 (1-xi)*(1-2*xi) xi*(2*xi-1) 0 4*xi*(1-xi) 0];
        l=sqrt((x(3)-x(2))^2+(y(3)-y(2))^2);
        elseif edgeno==3
            N=[xi*(2*xi-1) 0 (1-xi)*(1-2*xi) 0 0 4*xi*(1-xi)];
            l=sqrt((x(1)-x(3))^2+(y(1)-y(3))^2);
end

fe=N'*qn*l/2;