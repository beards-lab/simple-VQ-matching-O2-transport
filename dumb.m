clear; close all; clc;

x = 1:10;

y = stepz(x);


figure;
plot(x,y,'o')
grid on








function y = stepz(x)
y = zeros(1,length(x));
for i = 1:length(x)
    for j = 1:x(i)
        y(i) = y(i) + j;
    end
end
end