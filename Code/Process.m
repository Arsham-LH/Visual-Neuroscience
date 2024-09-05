%% Q3: Behavioral Data VS Coherence
clear;clc;

coh_levels = [0.032, 0.064, 0.128, 0.256];
subj_arr = ["Subj1", "Subj2"];
block_arr = 1:8;

N = length(subj_arr)*length(block_arr);
mean_RT_mat = zeros(length(coh_levels), N);
acc_mat = zeros(length(coh_levels), N);

j = 0;
for subj = subj_arr
    for block = block_arr
        j = j+1;
        path = sprintf('RDM/%s_block_%d.mat', subj, block);
        if ~exist(path)
            continue;
        end
        load(path);
        result = data.result;
        for i = 1:length(coh_levels)
            slct_result = result(result(:,2)==coh_levels(i) & result(:,7)==1 , :);
            mean_RT_mat(i,j) = mean(slct_result(:,6));
            N = size(slct_result, 1);
            acc_mat(i,j) = sum(slct_result(:,5)) / N;
        end
    end
end

mean_RT = mean(mean_RT_mat, 2);
std_RT = std(mean_RT_mat,[], 2);
mean_acc = mean(acc_mat, 2);
std_acc = std(acc_mat, [], 2);

figure;
sgtitle("Behavioral Data");

subplot(121);
errorbar(coh_levels*100, mean_acc, std_acc ./ sqrt(N), ".", "LineWidth", 1, "Color", "black");
xlabel("Motion Strength (%Coh)");
ylabel("Probability Correct");
xticks(coh_levels*100);
set(gca, 'XScale', 'log');
xlim([1.6, 51.2]);

subplot(122);
errorbar(coh_levels*100, mean_RT*1000, std_RT*1000 ./ sqrt(N), ".", "LineWidth", 1, "Color", "black");
xlabel("Motion Strength (%Coh)");
ylabel("Reaction Time (ms)");
xticks(coh_levels*100);
set(gca, 'XScale', 'log');
xlim([1.6, 51.2]);


%% Q4: Extract Phases
clear;clc;

coh_levels = [0.032, 0.064, 0.128, 0.256];
subj_arr = ["Subj1", "Subj2"];
block_arr = 1:8;
tot_data_P1 = [];
tot_data_P2 = [];
tot_data_P3 = [];

for subj = subj_arr
    for block = block_arr
        path = sprintf('RDM/%s_block_%d.mat', subj, block);
        if ~exist(path)
            continue;
        end
        load(path);
        if block <= 2
            tot_data_P1 = [tot_data_P1;data];
        elseif block >= 3 && block <= 6
            tot_data_P2 = [tot_data_P2;data];
        elseif block >= 7
            tot_data_P3 = [tot_data_P3;data];
        end
    end
end

% Total number of datasets for each phase:
N_P1 = length(tot_data_P1);
N_P2 = length(tot_data_P2);
N_P3 = length(tot_data_P3);

%% Q5: Behavioral Data VS Phase

coh_levels = [0.032, 0.064, 0.128, 0.256];
subj_arr = ["Subj1", "Subj2"];
block_arr = 1:8;

mean_RT_mat_P1 = zeros(N_P1, 1);
mean_RT_mat_P2 = zeros(N_P2, 1);
mean_RT_mat_P3 = zeros(N_P3, 1);
acc_mat_P1 = zeros(N_P1, 1);
acc_mat_P2 = zeros(N_P2, 1);
acc_mat_P3 = zeros(N_P3, 1);

for i = 1:N_P1
    data = tot_data_P1(i);
    result = data.result;
    slct_result = result(result(:,7)==1 , :);
    mean_RT_mat_P1(i) = mean(slct_result(:,6));
    N = size(slct_result, 1);
    acc_mat_P1(i) = sum(slct_result(:,5)) / N;
end
for i = 1:N_P2
    data = tot_data_P2(i);
    result = data.result;
    slct_result = result(result(:,7)==1 , :);
    mean_RT_mat_P2(i) = mean(slct_result(:,6));
    N = size(slct_result, 1);
    acc_mat_P2(i) = sum(slct_result(:,5)) / N;
end
for i = 1:N_P3
    data = tot_data_P3(i);
    result = data.result;
    slct_result = result(result(:,7)==1 , :);
    mean_RT_mat_P3(i) = mean(slct_result(:,6));
    N = size(slct_result, 1);
    acc_mat_P3(i) = sum(slct_result(:,5)) / N;
end


mean_RT = [mean(mean_RT_mat_P1), mean(mean_RT_mat_P2), mean(mean_RT_mat_P3)];
std_RT = [std(mean_RT_mat_P1), std(mean_RT_mat_P2), std(mean_RT_mat_P3)];
mean_acc = [mean(acc_mat_P1), mean(acc_mat_P2), mean(acc_mat_P3)];
std_acc = [std(acc_mat_P1), std(acc_mat_P2), std(acc_mat_P3)];
N_arr = [N_P1, N_P2, N_P3];



figure;
sgtitle("Behavioral Data")
subplot(121);
bar(1:3, mean_acc);
hold on;
errorbar(1:3, mean_acc, std_acc ./ sqrt(N_arr), ".", "LineWidth", 1, "Color", "black");
title("Accuracy");
xlabel("Phase");
ylabel("Accuracy");

subplot(122);
bar(1:3, mean_RT);
hold on;
errorbar(1:3, mean_RT, std_RT ./ sqrt(N_arr), ".", "LineWidth", 1, "Color", "black");
title("Reaction Time");
xlabel("Phase");
ylabel("Reaction Time (s)");
ylim([0, 1]);



%************* Required variables for Q7:
tot_mean_RT = mean_RT;
tot_mean_acc = mean_acc;


%% Q6: Write required matrix for each phase
clc;

printMat_P1 = [];
printMat_P2 = [];
printMat_P3 = [];

for i = 1:N_P1
    result = tot_data_P1.result;
    slct_result = result(result(:,7)==1 , :);
    acc_RT_mat = slct_result(:, [5,6]);
    printMat_P1 = [printMat_P1;acc_RT_mat];
end
for i = 1:N_P2
    result = tot_data_P2.result;
    slct_result = result(result(:,7)==1 , :);
    acc_RT_mat = slct_result(:, [5,6]);
    printMat_P2 = [printMat_P2;acc_RT_mat];
end
for i = 1:N_P3
    result = tot_data_P3.result;
    slct_result = result(result(:,7)==1 , :);
    acc_RT_mat = slct_result(:, [5,6]);
    printMat_P3 = [printMat_P3;acc_RT_mat];
end

% writematrix(printMat_P1,'DDM/phase1.dat','Delimiter','tab');
% writematrix(printMat_P2,'DDM/phase2.dat','Delimiter','tab');
% writematrix(printMat_P3,'DDM/phase3.dat','Delimiter','tab');

%% Q6: Plot DDM Parameters
clc;

DDM_mat = load("DDM/exp1Edit.logt");

a_arr = DDM_mat(:, 2); % Decision Bound
v_arr = DDM_mat(:, 3); % Drift Rate
t0_arr = DDM_mat(:, 4); % Non-decision Time


figure;
sgtitle("DDM Parameters");
subplot(131);
plot(1:3, v_arr, "Linewidth", 1);
hold on;
plot(1:3, v_arr, ".", "MarkerSize", 10, "Linewidth", 1);
title("Drift Rate");
xlabel("Phase");
ylabel("Drift Rate");
xlim([0.5, 3.5]);
xticks(1:3);

subplot(132);
plot(1:3, a_arr, "Linewidth", 1);
hold on;
plot(1:3, a_arr, ".", "MarkerSize", 10, "Linewidth", 1);
title("Decision Bound");
xlabel("Phase");
ylabel("Decision Bound");
xlim([0.5, 3.5]);
xticks(1:3);

subplot(133);
plot(1:3, t0_arr, "Linewidth", 1);
hold on;
plot(1:3, t0_arr, ".", "MarkerSize", 10, "Linewidth", 1);
title("Non Decision Time");
xlabel("Phase");
ylabel("Non Dec. time");
xlim([0.5, 3.5]);
xticks(1:3);

%% Q7: Fit Wang Model
clc;

% Choose the phase you want to fit model on (1 or 3):
slct_phase = 3;

behav_struct.tot_mean_RT = tot_mean_RT;
behav_struct.tot_mean_acc = tot_mean_acc;
behav_struct.printMat_P1 = printMat_P1;
behav_struct.printMat_P2 = printMat_P2;
behav_struct.printMat_P3 = printMat_P3;

CalCostFixedBehav = @(freeParams) CalCost(freeParams, behav_struct, slct_phase);

%Initial values for u0, thr:
freeParams0 = [17, 0.13];

options = optimset('MaxIter', 50);
[freeParams, min_cost] = fminsearch(CalCostFixedBehav, freeParams0, options);

fprintf('Phase %d: u0 = %d & thr = %d --> min_cost: %d', slct_phase, round(freeParams(1),2), round(freeParams(2),2), round(min_cost,2));

%% Q7: Plot Changes in Wang Model Parameters

u0_arr = [10.5, 13.4591];
thr_arr = [0.1687, 0.1829];

figure;
bar([1,2], [diff(u0_arr), diff(thr_arr)], "Linewidth", 1);
title("Wang Model Parameters");
xlabel("Parameters");
ylabel("Difference (Phase3 - Phase1)");
xticklabels(["u0","Threshold"]);


