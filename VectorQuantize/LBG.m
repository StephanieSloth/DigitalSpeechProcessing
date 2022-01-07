%% Vector Quantization
% 1.preparation:training set
% 2.clustering algorithm:L training data to M-element codebook
% 3.Nearest-neighbor codebook search

%% 1. preparation:training set
clc;clear;
% training set includes 2-Dimensional data of L(128)
trainset = load('training.dat');
L = length(trainset);
M = 4; % codebook size

%% 2. clustering algorithm:L training data to M-element codebook
% a. codebook initialization
center = sum(trainset)/L; % 1-vector codebook
e = 0.03; % 0.01<epsilon<0.05

for m=1:log2(M) %split log2(M) time
    D0 = 1; %D'
    D1 = 0; %D
    center = [center*(1+e);center*(1-e)];
    while abs(D1-D0)/D0>e % converge
        D0 = D1;
        % Nearest-Neighbor Search
        dmin = zeros(L,1);
        for i = 1:L % all train set
            for j = 1:length(center)
                % distance/distortion measure)
                d(i,j)=norm(trainset(i,:)-center(j,:));
                dmin_set = find(d(i,:)==min(d(i,:)));
                dmin(i) =dmin_set(1);
            end
        end

        % classify and center update
        for j = 1:length(center)
            cluster= find(dmin==j);
            cluster_set = trainset(cluster,:);
            center(j,:) = sum(cluster_set)/length(cluster_set);
        end

        % compute distortion
        distort = 0;
        for j = 1:length(center)
            cluster= find(dmin==j);
            cluster_set = trainset(cluster,:);
            dis = 0;
            for k = 1:length(cluster_set)
                dis = dis+norm(cluster_set(k,:)-center(j,:));
            end
            distort = distort+dis;
        end
        D1 = distort/L/2^m;
    end
end

%% 3. Nearest-neighbor codebook search
testset = load('to_be_quantized.dat');
L_test = length(testset);
dmin_test = zeros(L_test,1);
for i = 1:L_test % all train set
    for j = 1:length(center)
        % distance/distortion measure)
        dt(i,j)=norm(trainset(i,:)-center(j,:));
        dmin_test_set = find(dt(i,:)==min(dt(i,:)));
        dmin_test(i) =dmin_test_set(1);
    end
end
