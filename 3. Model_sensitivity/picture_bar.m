function picture_bar()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Plot: Sensitivity Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


picture_width = 26;
picture_hight = 14;
word = 14;
Xlabel = {'$$\theta_{D_0}$$';'$$\theta_{T_{N4}}$$';'$$\theta_{T_{N8}}$$';'$$\theta_M$$';
    '$$\kappa_{T_h}$$';'$$\kappa_{T_r}$$';'$$\kappa_{T_c}$$';'$$\beta_{T_h}$$';'$$\beta_{T_r}$$';
    '$$\beta_{T_c}$$'};
Ylabel = {'$$T_h$$' '$$T_r$$' '$$T_c$$' '$$TAM$$'};

%% load data

file=strcat('./data_prcc/','model_prcc.mat');
A = cell2mat(struct2cell(load(file)));
Color(1,:) = [152 251 152]/255;
Color(2,:) = [186 85 211]/255;
Color(3,:) = [255 99 71]/255;
Color(4,:) = [79 148 205]/255;
[~,c] = size(A');
figure(1)
for i =1:4
    c1=strcat(Ylabel(i)); 
    subplot(2,2,i)
    B = abs(A(:,i));
    bar(B,'FaceColor',Color(i,:));
    ax = gca;
    ax.XTick=1:1:c;
    ax.XLim=[0 c+0.5];
    ax.YLim=[1e-7,10];
    ax.TickLabelInterpreter='latex';
    ax.XTickLabel=Xlabel;
    ax.YScale = 'log';
    ax.Box='on';
    ax.FontName='Arial';
    ax.FontWeight='bold';
    set(gca,'FontSize',word);
    legend(c1, 'Location', 'northeast','Interpreter', 'latex');
end
    set(gcf,'unit','centimeters','position',[0.5 0.5 picture_width picture_hight]);
    file = strcat('picture_prcc/prcc.png');
    print('-dpng','-r600',file)

end