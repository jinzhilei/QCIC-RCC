function picture_basic()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot: cells and cytokines curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set: parameter
global par
parameter_basic();

%% cell
figure(1)
file = ['./basic_data/','cell.mat'];

% plot cell curve
A=cell2mat(struct2cell(load(file)));
t = 1:1:366;

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
end

hold on

% plot TN4 and TN8 in peripheral blood
S = [8,17];
for i = 1 : 2
    subplot(5,6,S(i))
    x_value = 0; 
    y_value = 3e8; 
    line([x_value, x_value+366], [y_value, y_value], 'Color', 'r', 'LineStyle', '--'); 
end

% plot K_D
subplot(5,6,4)
x_value = 0; 
y_value = par.K_D; 
line([x_value, x_value+366], [y_value, y_value], 'Color', 'r', 'LineStyle', '--');

% Set: image output
set(gcf,'unit','centimeters','position',[2 2 30 15])
print('-dpng','-r600','basic_picture\cell.png')


%%  chemokine
figure(2)
file=['./basic_data/','chemokine.mat'];

% plot chemokine curve
B=cell2mat(struct2cell(load(file)));
t = 1:1:366;
Y = {'I2_C','I10_C','I12_C','Igamma_C','Tbeta_C','L19_C',...
    'L21_C','X_D','C9_D','C10_D','L2_D','L17_D','L20_D','L22_D','TA'};
Z = [par.K_I2 par.K_I10 par.K_I12 par.K_Igamma par.K_Tbeta par.K_L19 par.K_L21 par.X par.K_C9 par.K_C10 par.K_L2 par.K_L17 par.K_L20 par.K_L22 par.K_TA];
for i = 1:1:15
    subplot(3,5,i);
    y = B(:,i);
    plot(t,y);
    ylabel(Y(i));
    hold on
    x_value = 0; 
    y_value = Z(i); 
    line([x_value, x_value+366], [y_value, y_value], 'Color', 'r', 'LineStyle', '--');
    hold off;
end

% Set: image output
set(gcf,'unit','centimeters','position',[2 2 35 15])
print('-dpng','-r600','basic_picture\chemokine.png')

end