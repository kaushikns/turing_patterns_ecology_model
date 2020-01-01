clearvars;
clc;
A = 3;
d = 0.3;
g = 0.8;
n = 100;
dl = 0.66/n;
l(1) = 0.14;
for i = 1:n
    r1(i) = (2*A*l(i) + (d+g)*(2*d+g) -(2*d+g)*sqrt(4*A*l(i)+(d+g)^2))/(2*l(i));
    l(i+1) = l(i) + dl;
end
l(i) = [];
area(l,r1, 'FaceColor',[.8 1 1]);
xlabel('\lambda');
ylabel('r');
xlim([0 0.8]);
ylim([0 1.0803]);

