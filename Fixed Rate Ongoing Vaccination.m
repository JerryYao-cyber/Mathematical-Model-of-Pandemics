% fixed number of vaccination for each day

clear all
close all 
% 7 days to recover 
% reproduction number 2
b = 1/7; 
a = 2*b;

% population
N = 3000; 
I = 3;
R = 0;
V = 0;
S = N-I-R-V;

% unit of time: day 
tmax = 200; % days
dt = 1/24; % 1 hour
clockmax = tmax/dt; 

%vaccination start day 0 to 0.02 
Vs = 50; 
Vpmax = 0.12;
dvp = 0.001;
stepmax = Vpmax/dvp;

Vrsave = zeros(1,stepmax+1); % vaccination proportion each day 
tIsave = zeros(1,stepmax+1); % total infection 
rIpsave = zeros(1,stepmax+1); % reduced infection proportion comparing to no vaccination 
%drIsave = zeros(1, stepmax); % difference of dI/ dvr
tpeak = zeros(1,stepmax+1); % time to reach the peak 

tsave = zeros(1,clockmax+1); 
Ssave = zeros(stepmax+1,clockmax+1);
Isave = zeros(stepmax+1,clockmax+1);
Rsave = zeros(stepmax+1,clockmax+1);
Vsave = zeros(stepmax+1,clockmax+1);
Vnewsave = zeros(stepmax+1,clockmax+1);


%Vr = 0, 0.001, 0.002,...,0.02
for Vr = 0.000:dvp:Vpmax % as vaccination rate varies

    column = int8((Vr/dvp)+1);
    Vrsave(column) = Vr; %  vaccination rate (as x-axis)

    tsave(column,1) = 0;
    Ssave(column,1) = S;
    Isave(column,1) = I;
    Rsave(column,1) = R;
    Vsave(column,1) = 0;

    for i = 1:(clockmax)
        if i*dt < Vs
            Vnewsave(column,i) = 0;
        else 
            if mod(i*dt,1)==0
                Vnewsave(column,i) = Vr * N;
            else
                Vnewsave(column,i) = 0;
            end
        end
    end

    for clock = 1: clockmax
        tsave(clock+1) = clock * dt;
        if clock > 1
            if Vnewsave(column,clock) > Ssave(column,clock-1)
                Vnewsave(column,clock) = Ssave(column,clock-1);
            end
        end
        I_new = Ssave(column,clock) * a * (Isave(column,clock)/N) * dt;
        R_new = Isave(column,clock) * b * dt;
        Ssave(column,clock+1) = Ssave(column,clock) - I_new - Vnewsave(column,clock);
        Isave(column,clock+1) = Isave(column,clock) + I_new - R_new;
        Rsave(column,clock+1) = Rsave(column,clock) + R_new;
        Vsave(column,clock+1) = Vsave(column,clock) + Vnewsave(column,clock);
    end
                                            
tIsave(column) = max(Rsave(column,:)); % the total infected people
rIpsave(column) = (tIsave(1) - tIsave(column))/N; %reduced infection proportion
end


Imaxsave = zeros(1,stepmax+1);
clockImaxsave = zeros(1,stepmax+1);

for step = 1:(stepmax+1)
    [Imaxsave(step), clockImaxsave(step)] = max(Isave(step,:));
end

subplot(4,1,1)
plot(Vrsave,tIsave)
title('Total Infected People')

subplot(4,1,2)
plot(Vrsave,rIpsave)
title('Reduced Infection Proportion')

subplot(4,1,3)
plot(Vrsave,Imaxsave)
title('Infection Peak')

subplot(4,1,4)
plot(Vrsave,clockImaxsave*dt)
title('Time when Infection Reaches the Peak')