function Data_processing()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% This function is used to process real patient data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load immunogenomic data

filePath  = 'Immunogenomic.xlsx';
True_data  = table2array(readtable(filePath));
CD8CD4_data=True_data(:,1);
CD4Treg_data=True_data(:,2);
TregTAM_data=True_data(:,3);

%% Save true data

file='./VP_immune/True_data_log.mat';
save(file,'True_data');

%% Data ksdensity

[PDF_CD8CD4,vector_data_CD8CD4]=ksdensity(CD8CD4_data);
[PDF_CD4Treg,vector_data_CD4Treg]=ksdensity(CD4Treg_data);
[PDF_TregTAM,vector_data_TregTAM]=ksdensity(TregTAM_data);

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
end
