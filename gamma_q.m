function fe=gamma_q(xi,coord,qn,edgeno)
x=coord(:,1);
y=coord(:,2);

if edgeno == 1
    l = sqrt((x(2)-x(1))^2 + (y(2)-y(1))^2);
    N1 = (1-xi)/2; N2 = (1+xi)/2; N3 = 0; N4 = 0;
elseif edgeno == 2
    l = sqrt((x(2)-x(3))^2 + (y(2)-y(3))^2);
    N1 = 0; N2 = (1-xi)/2; N3 = (1+xi)/2; N4 = 0;
elseif edgeno == 3
    l = sqrt((x(4)-x(3))^2 + (y(4)-y(3))^2);
    N1 =0; N2=0; N3 = (1+xi)/2; N4 = (1-xi)/2;
elseif edgeno == 4
    l = sqrt((x(4)-x(1))^2 + (y(4)-y(1))^2);
    N1 = (1-xi)/2; N2 = 0; N3 = 0; N4 = (1+xi)/2;
end

N = [N1, N2, N3,N4];

fe=N'*qn*l/2;