clearvars
close all

global values;

% 2000 simulation with SARS COV 2 parameters
impI = csvread("Out1.csv");
% impI2 = csvread("Out1.5.csv");
% impI3 = csvread("Out2.csv");
% impI4 = csvread("Out2.5.csv");impI5 = csvread("Out3.csv");
% impI6 = csvread("Out3.5.csv");impI7 = csvread("Out4.csv");
% impI8 = csvread("Out4.5.csv");impI9 = csvread("Out5.csv");


impI=impI';
impI2 = impI(1:end-1,:);

Mean = mean(impI2(end,2:end))
standardD = std(impI2(end,2:end));
[a,b]= max(impI2(end,2:end));

Data = impI2(end,1:end);
Data = Data';
% Data2 = Data (Data >= 0);

phatNBino = nbinfit(Data);
% phatpoiss = poissfit(impI2(end,2:end));

figure
h = histogram(impI2(end,1:end),78);
values = h.Values;
% legend('background','right ear-green','left ear-green')
title('Number of infected cells in 2000 of simulations')
ylabel('Number of simulation')
xlabel('Number of infected cells')

x3 = linspace(0,1);
figure

plot(x3,Gz(x3),'linewidth', 2);
hold on;
plot(x3,x3,'linewidth', 2);
ylabel('f(x)');
xlabel('x');
legend({'p.g.f.'},'Location','northwest');

x1 = 0;
x2 = Gz(x1);
format long
while abs(x2-x1) > 1e-10
    x1 = x2;
    x2 = Gz(x1);
end
x1
x2
% save('remove1.mat',x2)
%% Probability generating function
function out = Gz(x)
global values;

ps = values';
ps = ps / sum(ps);
[m,~] = size(ps);
% Gz = [1];
out = 0;
    for i = 1:m   
        out = out + ps(i) * x.^ (i-1);
    %    Gz = conv(Gz , ps(i)); 
       %conv as a polynomial multiplication
    end
end


