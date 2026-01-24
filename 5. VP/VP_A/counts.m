function  counts()
%% Set parameters

sample = 100;
z = 1000;

%% Load data

file_R='./VP_A_bar/VP_tumor_initialzation.mat';
r = cell2mat(struct2cell(load(file_R)));
file='./VP_A_bar/VPsample.mat';
you = cell2mat(struct2cell(load(file)));

%% Simulation calculation

L = length(you);
    CPSP = zeros(z,4);
    tumor = zeros(sample,1);
    ratio = zeros(sample,1);

% delete (gcp('nocreate'))            % Close the current parallel computing pool
% parpool(7);                         % Create a parallel pool

for j = 1: z
    disp(j)

    Sampling_Number = randi([1 L], sample, 1);
    tumor_0 = ((r(Sampling_Number).*2.572e-9.*0.8*3)./(0.37*4*3.14)).^(1/3);
    yout = you(Sampling_Number,:);
    for i = 1:sample
        tumor(i) = yout(i,25)+yout(i,26)+yout(i,27)+yout(i,28);
        tumor(i) = ((3*tumor(i).*2.572e-9.*0.8)./(0.37*4*3.14)).^(1/3);
        ratio(i) = (tumor(i)- tumor_0(i))./tumor_0(i)*100;

        if ratio(i) >= 20
            CPSP(j,1) = CPSP(j,1)+1;
        elseif -30 <= ratio(i) && ratio(i) < 20
            CPSP(j,2) = CPSP(j,2)+1;
        elseif -80 <= ratio(i) && ratio(i) <-30
            CPSP(j,3) = CPSP(j,3)+1;
        elseif -80 > ratio(i)
            CPSP(j,4) = CPSP(j,4)+1;
        end

    end
end

%% Save Result

CPSP = CPSP./sample.*100;
file='./VP_A_bar/CPSP';
save(file,'CPSP');

end

