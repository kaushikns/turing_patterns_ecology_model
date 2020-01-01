function [ err ] = lambda_c(l)
A = 3;
d = 0.3;
r = 0.5;
g = 0.8;
ds = 0.02;
di = 0.0005;
R = l*A/(d*(d+g));
H = l*r/(d*(d+g));
i = (d/(2*l))*(R - 1 - H + sqrt((R-1-H)^2 - 4*H));
s = A/(d + l*i);
j11 = -d - l*i;
j12 = -l*s;
j21 = l*i;
j22 = l*s - d - g;
lhs = ((j11*di + j22*ds)^2)/(4*di*ds);
rhs = j11*j22 - j12*j21;
err = abs(lhs-rhs);
end

