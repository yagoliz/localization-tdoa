
clear all
close all
clc

X_data_Example = [1+rand(1,50) 3+rand(1,50)]';
Y_data_Example = [1+rand(1,50) 4+rand(1,50)]';

epsilon = 0.7;
num_outliers_percent = 0.07;
num_mic_C_percent = 0.07;
Wc = 1;
W1 = 1;
W2 = 2;

subplot(2,1,1)
plot(X_data_Example,Y_data_Example,"*")
xlim([0 6])
ylim([0 6])
grid on
set(gca,'fontsize',24);

matrix_Example = [X_data_Example, Y_data_Example];
total_points = numel(X_data_Example);
num_outliers_per_cluste = max([1 floor(num_outliers_percent*total_points)]);
idx = dbscan(matrix_Example,epsilon,num_outliers_per_cluste);
Num_estimated_claster = max(idx);

C = zeros(1,Num_estimated_claster);
Var_1 = zeros(1,Num_estimated_claster);
Var_2 = zeros(1,Num_estimated_claster);
Mean_1 = zeros(1,Num_estimated_claster);
Mean_2 = zeros(1,Num_estimated_claster);

for ii = 1:Num_estimated_claster
    C(ii) = sum(idx == ii)./total_points;
    Var_1(ii) = var(Y_data_Example(idx==ii));
    Var_2(ii) = var(X_data_Example(idx==ii));
    Mean_1(ii) = mean(Y_data_Example(idx==ii));
    Mean_2(ii) = mean(X_data_Example(idx==ii));
end

if Num_estimated_claster == 1
    Choise_1 = Mean_1;
    Choise_2 = Mean_2;
else
    Likelihood = (Wc*C - W1*Var_1 - W2*Var_2);
    Likelihood(C<num_mic_C_percent) = -inf;
    Choise_1 = Mean_1(Likelihood == (max(Likelihood)));
    Choise_2 = Mean_2(Likelihood == (max(Likelihood)));
end

subplot(2,1,2)
for ii = 1:Num_estimated_claster
    plot(X_data_Example(ii == idx),Y_data_Example(ii == idx),"*")
    hold on
end
legend
xlim([0 6])
ylim([0 6])
grid on
set(gca,'fontsize',24);

