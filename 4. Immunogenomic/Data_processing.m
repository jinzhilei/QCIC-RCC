function Data_processing(sample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This function is used to process both virtual and real patient data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load immunogenomic data

filePath  = 'Immunogenomic.xlsx';
CellType  = table2array(readtable(filePath));
CD8CD4_data=CellType(:,1);
CD4Treg_data=CellType(:,2);
TregM_data=CellType(:,3);


%% Save true data

True_data  = zeros(959,3);
True_data(:,1) = CD8CD4_data';
True_data(:,2) = CD4Treg_data';
True_data(:,3) = TregM_data';

file='./VP_immune/True_data_log.mat';
save(file,'True_data');

%% Data ksdensity

[PDF_CD8CD4,vector_data_CD8CD4]=ksdensity(CD8CD4_data);
[PDF_CD4Treg,vector_data_CD4Treg]=ksdensity(CD4Treg_data);
[PDF_TregTAM,vector_data_TregTAM]=ksdensity(TregM_data);

%% Save distribution

PDF=zeros(6,100);
PDF(1,:)=vector_data_CD8CD4;
PDF(2,:)=PDF_CD8CD4;
PDF(3,:)=vector_data_CD4Treg;
PDF(4,:)=PDF_CD4Treg;
PDF(5,:)=vector_data_TregTAM;
PDF(6,:)=PDF_TregTAM;

file='./VP_immune/True_data_ksdensity.mat';
save(file,'PDF');

%% Load Plausible patients data

file  = './VP_immune/VP_sample_immune.mat';
A     = cell2mat(struct2cell(load(file)));
CD4   = A(:,12) ;
CD8   = A(:,21);
Treg  = A(:,15);
TAM   = A(:,24);

B = zeros(sample,3);
B(:,1) = CD8./CD4;
B(:,2) = CD4./Treg;
B(:,3) = Treg./TAM;
B = log(B);

file='./VP_immune/VP_data_log.mat';
save(file,'B');
