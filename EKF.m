% Surface Temperature Estimation of Li-Ion Battery using extended Kalman Filter
% Both State and Measurement have nonlinear models
% (State Model parameters based on Cycle 50 and Measurement Model
% parameters based on Cycle 200)
% Using B0005-Cycle 450 Data

%clear all;
clc
close all;

% Li-Ion Battery Data

%TempdegC = xlsread('B0005.xls','Cycle 450','G:G');
%CurrentA = xlsread('B0005.xls','Cycle 450','E:E');
%VoltageA = xlsread('B0005.xls','Cycle 450','C:C');
T_prev = TempdegC(3:268,1);
T = TempdegC(4:269,1);
Current = CurrentA(3:268,1);
V = VoltageV(3:269,1);

n=266;

%Converting Column vector to Row vector
T = T';
Current = Current';
V = V';


q=0.3; %Just assuming
r=0.01;


% Do the initialization
%Our x_truth is T(1:n)
x_truth(1) = T(1);
z(1) = V(1); %Our actual measurement is the Voltage array  

%Initialization
x_bar(1) = -(log(z(1))-30.6361)/0.82726;
x_hat(1) = x_bar(1);

%Initializing Variances

Pbar(1) = (x_truth(1)-x_bar(1))^2;
Phat(1) = Pbar(1);

% Simulate the Truth

for i = 2:n

    %%% 2. Simulate system dynamics by adding noise

    x_truth(i) = T(i);

    %%% 3. Simulate the measurement by adding noise
    z(i) = V(i);  

end



% Do the prediction

for i = 2:n
    
    x_bar(i) =  0.7993*x_hat(i-1)+0.05315*Current(i)^2+0.00312*x_hat(i-1)^2-0.0003826*x_hat(i-1)*Current(i)^2+3.108;%State equation evaluated at previous estimate x_hat(i-1)
    
    zbar(i) = -8.947e-15*exp(0.8136*x_bar(i))+5.534*exp(-0.01366*x_bar(i)); % the measurement equation- z(k)=h(k,x) evaluated at x_bar(i)-state prediction
    
    z_tilda = z(i) - zbar(i);
    %Calculating Gradients
    fG = 0.00624*x_hat(i-1)-0.0003826*Current(i)^2+0.7993;  %Gradient of f evaluated at previous estimate x_hat(i-1)
    hG = -7.2793e-15*exp(0.8136*x_bar(i)-0.075594*exp(-0.01366*x_bar(i))); %Gradient of h evaluated at x_bar(i)
    
    %State Prediction Covariance
    Pbar(i) = fG * Phat(i-1) * fG' + q;
    
    %Innovation Covariance or
    S = hG * Pbar(i) * hG' + r;
    
    W = Pbar(i) * hG' * inv(S);
    
    x_hat(i) = x_bar(i) + W * z_tilda;
    
    Phat(i) = Pbar(i) - W * S * W';
end

   
i = 1:n;


MSE = (x_truth-x_hat).^2/n;

i = 1:n;


figure

plot(i, x_truth(:),'b-', i , x_hat(:),'r-')
xlabel('time');
ylabel('Actual versus Estimated Temperature')
ylim([0 45])
title('EKF-Model 1')

figure

plot(i, ((x_truth(:)-x_hat(:)).^2)/n,'b-')
xlabel('time');
ylabel('MSE')
ylim([0 0.6])
title('Mean Square Error')




    

