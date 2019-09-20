function [ output_args ] = davies_bouldin( input_args )
%DAVIES_BOULDIN Summary of this function goes here
%   Detailed explanation goes here
    parfor j=1:columns
        %%Davies–Bouldin index
        MDB=zeros(3);
        R=zeros(3);
        D=zeros(3,1);
        %calculating avg of each cluster
        avg1(j) = mean(M1(:,j));
        avg0(j) = mean(M0(:,j));
        avg_1(j) = mean(M_1(:,j));

        %calculating intracluster distance
        S1=sqrt(sum(power(M1(:,j)-avg1(j),2))/n1);
        S0=sqrt(sum(power(M0(:,j)-avg0(j),2))/n0);
        S_1=sqrt(sum(power(M_1(:,j)-avg_1(j),2))/n_1);

        %calculating intercluster distance
        MDB(1,2)=abs(avg1(j)-avg0(j));
        MDB(2,1)=MDB(1,2);
        MDB(1,3)=abs(avg1(j)-avg_1(j));
        MDB(3,1)=MDB(1,3);
        MDB(2,3)=abs(avg0(j)-avg_1(j));
        MDB(3,2)=MDB(2,3);

        %calculating intercluster distance
        R(1,2)=(S1+S0)/MDB(1,2);
        R(2,1)=R(1,2);
        R(1,3)=(S1+S_1)/MDB(1,3);
        R(3,1)=R(1,3);
        R(2,3)=(S0+S_1)/MDB(2,3);
        R(3,2)=R(2,3);

        %saving current gene's index
        D(1)=max(R(1,2),R(1,3));
        D(2)=max(R(2,1),R(2,3));
        D(3)=max(R(3,1),R(3,2));
        DB(j)=sum(D)/3;
    %     sumMDB(j)=sum(sum(MDB));
    %     if ((sumMDB(j)>200)&&(DB(j)<50))||((sumMDB(j)>100)&&(DB(j)<8))
    %     if DB(j)<3.955
    %         savedGenes(j)=1;
    %     end

        %%
    %       s=silhouette(M_valid(:,j), r(rows,:));
    %       for i=1:203
    %           if s(i)<-1||s(i)>1
    %               s(i)=0;
    %           end
    %       end
    %       avgSilhouette(j)=sum(s);sum(M_1(:,j)-avg_1(j))/n_1
    %       
    end
end

