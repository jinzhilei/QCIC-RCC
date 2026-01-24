function counts_recurrence_rate()
sample =100;
Time = 360;
A = zeros(3,1);
B = zeros(3,1);
tumor_0 = zeros(sample,1);
for b  =1:3
    for i =1:sample
        file=['./Adaptive_therapy/Adaptive_therapy_7/',num2str(b),'/VPsample_AB',num2str(i),'.mat'];
        cell=cell2mat(struct2cell(load(file)));
        tumor_0(i) = ((3*cell(1,25).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
        tumor(:,i) = cell(:,25)+cell(:,26)+cell(:,27)+cell(:,28);
        tumor(:,i) = ((3*tumor(:,i).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
        ratio(:,i) = tumor(:,i)- tumor_0(i);
        if ratio(Time,i) >= 0.2
            A(b,1) = A(b,1)+1;
        end
    end
end

A = A./sample;

for b  =1:3
    for i =1:sample
        file=['./Adaptive_therapy/Adaptive_therapy_21/',num2str(b),'/VPsample_AB',num2str(i),'.mat'];
        cell=cell2mat(struct2cell(load(file)));
        tumor_0(i) = ((3*cell(1,25).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
        tumor(:,i) = cell(:,25)+cell(:,26)+cell(:,27)+cell(:,28);
        tumor(:,i) = ((3*tumor(:,i).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
        ratio(:,i) = tumor(:,i)- tumor_0(i);
        if ratio(Time,i) >= 0
            B(b,1) = B(b,1)+1;
        end
    end
end

B = B./sample;

file=['Recurrence_rate/treat_7/','Recurrence_rate.mat'];
save(file,'A');

file=['Recurrence_rate/treat_21/','Recurrence_rate.mat'];
save(file,'B');
end

