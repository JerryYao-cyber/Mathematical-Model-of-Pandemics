
clear all
close all
% 7 days to recover
% reproduction number 2.3
% unit of time is one day

b = 1/7;
a = 2.3*b;

N = 10000;
I = 3;
R = 0;
% S =? 

tmax = 500; % days
dt = 20; % 1 hour
clockmax = tmax/dt; 

Vp_max = 0.8; % 80percent of the population is inherently immune
dvp = 0.05; % change in percentage
stepmax = Vp_max/dvp;

Vrsave = zeros(1,stepmax+1); 
tIsave = zeros(1,stepmax+1); % total number of infected ppl which is R max
tIpsave = zeros(1,stepmax+1); % total proportion of the S(t=0) that is infected when immune propotion vary
tpeaksave = zeros(1,stepmax+1); % time takes to reach the peak 
unIsave = zeros(1,stepmax+1);

tsave = zeros(1,clockmax+1);  
Isave = zeros(stepmax+1,clockmax+1);
Rsave = zeros(stepmax+1,clockmax+1);
Ssave = zeros(stepmax+1,clockmax+1);


for step = 1:(stepmax+1)
    V = (step-1) * dvp * N;
    Ssave(step,1) = N - I - R - V;
    Isave(step,1) = I;
    Rsave(step,1) = R;

    tsave(1)= 0;
    for clock = 1: clockmax
        tsave(clock+1) = clock * dt;
        I_new = Ssave(step,clock) * a * (Isave(step,clock)/N) * dt;
        R_new = Isave(step,clock) * b * dt;
        Ssave(step,clock+1) = Ssave(step,clock) - I_new;
        Isave(step,clock+1) = Isave(step,clock) + I_new - R_new;
        Rsave(step,clock+1) = Rsave(step,clock)+ R_new;
    end
    Vrsave(step) = (step-1) * dvp;
    tIsave(step) = max(Rsave(step,:));
    tIpsave(step) = tIsave(step)/Ssave(step,1); % total proportion of the S(t=0) that is infected when immune propotion vary
    unIsave(step) = N - V - max(Rsave(step,:));%not infected not immune
end

Imaxsave = zeros(1,9);
clockImaxsave = zeros(1,9);
for index = 1:9 
    [Imaxsave(index),clockImaxsave(index)]= max(Isave(index,:));
end

%tIsave = max(Rsave.');
%tIpsave = tIsave/(Ssave(:,1).');

subplot(4,1,1)
plot(Vrsave, tIsave, 'r')
%xlabel('Innate Immune Percentage')
%ylabel('Total Infection')
title('Total # of infected people at the end of the pandemics vs. Initial immune percentage')

subplot(4,1,2)
plot(Vrsave, tIpsave, 'b')
%xlabel('Innate Immune Percentage')
%ylabel('Total Infection Percentage')
title('Total fraction # of infected people at the end of the pandemics vs. Initial immune percentage')

subplot(4,1,3)
plot(Vrsave, clockImaxsave*dt, 'g')
%xlabel('Innate Immune Percentage')
%ylabel('Time of Infection Peak')
title('# of days for infected population to reach its maximum vs. Initial immune percentage')

subplot(4,1,4)
plot(Vrsave, unIsave,'m' )
%xlabel('Innate Immune Percentage')
%ylabel('Uninfected S')
title('# of uninfected susceptible people vs. Initial immune percentage')


% exceed N



