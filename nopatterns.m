clearvars;
clc;

% No Pattern formation
A = 3;
d = 0.3;
r = 0.5;
g = 0.8;
ds = 0.02;
di = 0.0005;
l = 0.65;
dx = 0.05;
dy = 0.05;
dt = 0.01;
R = l*A/(d*(d+g));
H = l*r/(d*(d+g));
i0 = (d/(2*l))*(R - 1 - H + sqrt((R-1-H)^2 - 4*H));
s0 = A/(d + l*i0);
nx = 5/dx;
ny = 5/dy;
nt = 20/dt;
i(nt,nx,ny) = 0;
s(nt,nx,ny) = 0;
for p=1:nx
    for q=1:ny
        i(1,p,q) = i0*(1+(-1)^p*(-1)^q*rand);
        s(1,p,q) = s0*(1+(-1)^p*(-1)^q*rand);
    end
end
for m=1:nt-1
    for p=2:nx-1
        for q=2:ny-1
            i(m+1,p,q) = i(m,p,q) + dt*(l*s(m,p,q)*i(m,p,q) - (d+g)*i(m,p,q) - r + di*((i(m,p,q+1)+i(m,p,q-1)-2*i(m,p,q))/(dy*dy) + (i(m,p+1,q)+i(m,p-1,q)-2*i(m,p,q))/(dx*dx)));
            s(m+1,p,q) = s(m,p,q) + dt*(-l*s(m,p,q)*i(m,p,q) - d*s(m,p,q) + A + ds*((s(m,p,q+1)+s(m,p,q-1)-2*s(m,p,q))/(dy*dy) + (s(m,p+1,q)+s(m,p-1,q)-2*s(m,p,q))/(dx*dx)));
        end
        i(m+1,p,1) = i(m+1,p,2);
        i(m+1,p,ny) = i(m+1,p,ny-1);
        s(m+1,p,1) = s(m+1,p,2);
        s(m+1,p,ny) = s(m+1,p,ny-1);
    end
    for q=1:ny
        i(m+1,1,q) = i(m+1,2,q);
        i(m+1,nx,q) = i(m+1,nx-1,q);
        s(m+1,1,q) = s(m+1,2,q);
        s(m+1,nx,q) = s(m+1,nx-1,q); 
    end
    disp(m);
end
iplot(nx,ny) = 0;
for p=1:nx
    for q=1:ny
        iplot(p,q) = i(m,p,q);
    end
end
HeatMap(iplot);
