function cost = CalCost(freeParams, behav_struct, slct_phase)
    behav_mean_RT = behav_struct.tot_mean_RT(slct_phase);
    behav_mean_acc = behav_struct.tot_mean_acc(slct_phase);
    
    printMat = [];
    if slct_phase == 1
        printMat = behav_struct.printMat_P1;
    elseif slct_phase == 2
        printMat = behav_struct.printMat_P2;
    elseif slct_phase == 3
        printMat = behav_struct.printMat_P3;
    end

    C1 = WANG_E_edit(freeParams(1), freeParams(2));
    C1 = C1(:,[5,6]); %Keeping only RT & Acc
    model_mean_RT = mean(C1(:,1));
    model_mean_acc = mean(C1(:,2));
    
    RT_cost = (behav_mean_RT - model_mean_RT)^2;
    acc_cost = (behav_mean_acc - model_mean_acc)^2;

    %************ RT Histogram
    edges = 0:0.05:2;

    behav_hist = histogram(printMat(:,2), edges, 'Normalization', 'probability' );
    behav_vals = behav_hist.Values;
    model_hist = histogram(C1(:,1), edges, 'Normalization', 'probability' );
    model_vals = model_hist.Values;

    hist_cost = sum((behav_vals - model_vals).^2);
    cost = 200*acc_cost + 2*RT_cost + 15*hist_cost;
end





