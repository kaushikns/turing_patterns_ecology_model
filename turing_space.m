clearvars;
clc;

% Determination of lambda_c and kc
A = 3;
d = 0.3;
r = 0.5;
g = 0.8;
ds = 0.02;
di = 0.0005;
lc = 0.8;
lc = fsolve(@lambda_c,lc);
R = lc*A/(d*(d+g));
H = lc*r/(d*(d+g));
i = (d/(2*lc))*(R - 1 - H + sqrt((R-1-H)^2 - 4*H));
s = A/(d + lc*i);
j11 = -d - lc*i;
j12 = -lc*s;
j21 = lc*i;
j22 = lc*s - d - g;
kc = sqrt((j11*di + j22*ds)/(2*di*ds));

% Bifurcation diagram 1
n = 1;
l(n) = 0.4;
dl = 0.01;
while l(n) < 0.7
    R = l(n)*A/(d*(d+g));
    H = l(n)*r/(d*(d+g));
    i = (d/(2*l(n)))*(R - 1 - H + sqrt((R-1-H)^2 - 4*H));
    s = A/(d + l(n)*i);
    j11 = -d - l(n)*i;
    j12 = -l(n)*s;
    j21 = l(n)*i;
    j22 = l(n)*s - d - g;
    J = [j11 j12; j21 j22];
    k1(n) = sqrt((j11*di + j22*ds - sqrt((j11*di+j22*ds)^2 - 4*ds*di*det(J)))/(2*ds*di));
    k2(n) = sqrt((j11*di + j22*ds + sqrt((j11*di+j22*ds)^2 - 4*ds*di*det(J)))/(2*ds*di));
    n = n+1;
    l(n) = l(n-1) + dl;
end
l(n) = [];
plot(l,k1);
hold on;
plot(l,k2);

figure;

% Bifurcation diagram 2 & 3
dl = 0.1;
n = 1;
l1(n) = 0.3467;
dk_2 = 10;
while l1(n) < 0.9
    R = l1(n)*A/(d*(d+g));
    H = l1(n)*r/(d*(d+g));
    i = (d/(2*l1(n)))*(R - 1 - H + sqrt((R-1-H)^2 - 4*H));
    s = A/(d + l1(n)*i);
    j11 = -d - l1(n)*i;
    j12 = -l1(n)*s;
    j21 = l1(n)*i;
    j22 = l1(n)*s - d - g;
    J = [j11 j12; j21 j22];
    q = 1;
    k_2(q) = 0;
    while k_2(q) < 500
        h(n,q) = k_2(q)*k_2(q)*di*ds - k_2(q)*(j22*ds+j11*di) + det(J);
        f = [1, -(j11 - k_2(q)*(ds+di) + j22), h(n,q)];
        delta = roots(f);
        re_delta(n,q) = real(delta(2));
        q = q + 1;
        k_2(q) = k_2(q-1) + dk_2;
    end
    n = n+1;
    l1(n) = l1(n-1) + dl;
end
k_2(q) = [];
l1(n) = [];
q = 1;
k_2(q) = 0;
while k_2(q) < 500
    h(n,q) = 0;
    re_delta(n,q) = 0;
    q = q + 1;
    k_2(q) = k_2(q-1) + dk_2;
end
k_2(q) = [];
for p = 1:n-1
    plot(k_2,h(p,:));
    hold on;
end
p = p+1;
plot(k_2,h(p,:),'--');
figure;

for p=1:n-1
    plot(k_2,re_delta(p,:));
    hold on;
end
p = p +1;
plot(k_2,re_delta(p,:),'--');
