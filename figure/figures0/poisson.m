clf;
hold on 
x = 1:1:40;
lambda = 5;
a = 4;
x1 = 1:a:40;
y1 = poisspdf(x1,lambda*a)*a;
plot(1:size(x1,2),y1,':r*');

a = 2;
x1 = 1:a:40;
y1 = poisspdf(x1,lambda*a)*a;
plot(1:size(x1,2),y1,'g*');

a = 1;
x1 = 1:a:40;
y1 = poisspdf(x1,lambda*a)*a;
plot(1:size(x1,2),y1,'b');