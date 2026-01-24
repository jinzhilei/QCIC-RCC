function picture_basic()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot: cells and cytokines curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set: parameter

%% cell
figure(1)
for j = 1:100

    file=['./VP_A/','VPsample_A',num2str(j),'.mat'];

    % plot cell curve
    A=cell2mat(struct2cell(load(file)));
    t = 1:1:360;

    Y = {'D0_A','D0_B','D0_D','D_C','D_D','D_E',...
        'TN4_A','TN4_B','TN4_C','Th_B','Th_C','Th_D',...
        'Tr_B','Tr_C','Tr_D','TN8_A','TN8_B','TN8_C',...
        'Tc_B','Tc_C','Tc_D','M_A','M_B','TAM_D',...
        'X1_D','X2_D','X3_D','X4_D'};

    for i = 1:1:28
        subplot(5,6,i);
        y = A(:,i);
        plot(t,y);
        ylabel(Y(i));
        hold on
    end

    hold on

    % plot TN4 and TN8 in peripheral blood
    S = [8,17];
    for i = 1 : 2
        subplot(5,6,S(i))
        x_value = 0;
        y_value = 3e8;
        line([x_value, x_value+366], [y_value, y_value], 'Color', 'r', 'LineStyle', '--'); % 横线的长度和样式可以根据需要调整
    end
    hold on
    % plot K_D
    subplot(5,6,4)
    x_value = 0;
    line([x_value, x_value+366], [y_value, y_value], 'Color', 'r', 'LineStyle', '--'); 
    hold on
end

figure(2)
for j = 1:100

    file=['./VP_A/','VPsample_A',num2str(j),'.mat'];

    % plot cell curve
    A=cell2mat(struct2cell(load(file)));
    t = 1:1:360;
 plot(t,A(:,25));
 hold on
end

end