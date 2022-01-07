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
% draw 1
% figure(1);
% scatter(trainset(:,1),trainset(:,2));
% c0 = sum(trainset)/L; 
% hold on;
% scatter(c0(1),c0(2),200,'r','+');
% title('Centroid of Training Vectors');
%
e = 0.01; % 0.01<epsilon<0.05

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
       % draw 2
%     figure(2);
%    clu1= find(dmin==1);
%    clu1_set = trainset(clu1,:);
%    color1 = [0.8500 0.3250 0.0980];
%    scatter(clu1_set(:,1),clu1_set(:,2),[],color1);
%    hold on;
%    c1 = center(1,:);
%    scatter(c1(1),c1(2),200,'r','+');
%    clu2= find(dmin==2);
%    clu2_set = trainset(clu2,:);
%    color2 = [0 0.4470 0.7410];
%    scatter(clu2_set(:,1),clu2_set(:,2),[],color2,'s');
%    hold on;
%    c2 = center(2,:);
%    scatter(c2(1),c2(2),200,'r','+');
%    title('Two Vectors after Splitting');

 % draw 4
    
   clu1= find(dmin==1);
   clu1_set = trainset(clu1,:);
   color1 = [0.8500 0.3250 0.0980];
   c1 = center(1,:);

   clu2= find(dmin==2);
   clu2_set = trainset(clu2,:);
   color2 = [0.4660 0.6740 0.1880];
   c2 = center(2,:);

   clu3= find(dmin==3);
   clu3_set = trainset(clu3,:);
   color3 = [0 0.4470 0.7410];
   c3 = center(3,:);

   clu4= find(dmin==4);
   clu4_set = trainset(clu4,:);
   color4 = [0.4940 0.1840 0.5560];
   c4 = center(4,:);


   figure;
   scatter(clu1_set(:,1),clu1_set(:,2),[],color1,'o');
   hold on;
   scatter(clu2_set(:,1),clu2_set(:,2),[],color2,'s');
   hold on;
   scatter(clu3_set(:,1),clu3_set(:,2),[],color3,'d');
   hold on;
   scatter(clu4_set(:,1),clu4_set(:,2),[],color4,'^');
   scatter(c1(1),c1(2),200,'r','+');
   scatter(c2(1),c2(2),200,'r','+');
   scatter(c3(1),c3(2),200,'r','+');
   scatter(c4(1),c4(2),200,'r','+');
   title('Four Vectors after Splitting');

%% 3. Nearest-neighbor codebook search

testset = load('to_be_quantized.dat');
L_test = length(testset);
dmin_test = zeros(L_test,1);
for i = 1:L_test % all train set
    for j = 1:length(center)
        % distance/distortion measure)
        dt(i,j)=norm(testset(i,:)-center(j,:));
        dmin_test_set = find(dt(i,:)==min(dt(i,:)));
        dmin_test(i) =dmin_test_set(1);
    end
end

% compute distortion
        distort_test = 0;
        for j = 1:length(center)
            cluster= find(dmin_test==j);
            cluster_set = testset(cluster,:);
            dis = 0;
            for k = 1:length(cluster_set)
                dis = dis+norm(cluster_set(k,:)-center(j,:));
            end
            distort_test = distort_test+dis;
        end
        D1_test = distort_test/L/2^m;
% figure;
% scatter(testset(:,1),testset(:,2),'x');
    
   clu1= find(dmin_test==1);
   clu1_set = testset (clu1,:);
   color1 = [0.8500 0.3250 0.0980];
   c1 = center(1,:);

   clu2= find(dmin_test==2);
   clu2_set = testset (clu2,:);
   color2 = [0.4660 0.6740 0.1880];
   c2 = center(2,:);

   clu3= find(dmin_test==3);
   clu3_set = testset (clu3,:);
   color3 = [0 0.4470 0.7410];
   c3 = center(3,:);

   clu4= find(dmin_test==4);
   clu4_set =testset (clu4,:);
   color4 = [0.4940 0.1840 0.5560];
   c4 = center(4,:);


   figure;
   scatter(clu1_set(:,1),clu1_set(:,2),[],color1,'o');
   hold on;
   scatter(clu2_set(:,1),clu2_set(:,2),[],color2,'s');
   hold on;
   scatter(clu3_set(:,1),clu3_set(:,2),[],color3,'d');
   hold on;
   scatter(clu4_set(:,1),clu4_set(:,2),[],color4,'^');
   scatter(c1(1),c1(2),200,'r','+');
   scatter(c2(1),c2(2),200,'r','+');
   scatter(c3(1),c3(2),200,'r','+');
   scatter(c4(1),c4(2),200,'r','+');
   title('Test Data after Splitting');
