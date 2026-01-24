function  KL_JS_div()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% calculate: divergence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set parameter

score_KL = zeros(3,1);
score_JS = zeros(3,1);

%% Loda data
file = ['./VP_immune/','Number_VP_patients.mat'];
number=cell2mat(struct2cell(load(file)));
file=['./VP_immune/','VP_sample_immune.mat'];
A=cell2mat(struct2cell(load(file)));
CD4 = A(number,12) ;
CD8 = A(number,21);
Treg = A(number,15);
TAM = A(number,24);
B(:,1) = CD8./CD4;
B(:,2) = CD4./Treg;
B(:,3) = Treg./TAM;
B = log(B);

[PDF_VP(2,:),PDF_VP(1,:)]=ksdensity(B(:,1));
[PDF_VP(4,:),PDF_VP(3,:)]=ksdensity(B(:,2));
[PDF_VP(6,:),PDF_VP(5,:)]=ksdensity(B(:,3));

file = './VP_immune/True_data_ksdensity.mat';
PDF=cell2mat(struct2cell(load(file)));
PDF_P = PDF_VP;
for i = 2:2:6
    if any(PDF(i,:))
        PDF(i,:) = PDF(i,:)./sum(PDF(i,:));
    end
    if any(PDF_P(i,:))
        PDF_P(i,:) = PDF_P(i,:)./sum(PDF_P(i,:));
    end
    % Compute Kullback-Leibler Divergence
    score_KL(i/2) = sum(sum(PDF(i,:).* log(eps + PDF(i,:)./(PDF_P(i,:)+eps))));
    % Compute Jensen-Shannon Divergence
    score_JS(i/2) = (sum(sum(PDF(i,:).* log(eps + PDF(i,:)./((PDF(i,:)+PDF_P(i,:))./2+eps))))+sum(sum(PDF_P(i,:).* log(eps + PDF_P(i,:)./((PDF(i,:)+PDF_P(i,:))./2+eps)))))./2;
end

if PDF==PDF_VP
    score_KL=0;
    score_JS=0;
end

%% Save data
file='./VP_immune/score_KL.mat';
save(file,'score_KL');
file='./VP_immune/score_JS.mat';
save(file,'score_JS');
file='./VP_immune/PDF_VP.mat';
save(file,'PDF_VP');

end