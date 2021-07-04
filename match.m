% nm: the normal tissue sample matrix for patients
% tm: the tumor tissue sample matrix for patients

Disease_sig2=[];
mat1=nm-tm;
[~,ns]=size(mat1);
for i=1:12328
    temp=mat1(maps_LINCS_probes{i,3},:);
    if length(maps_LINCS_probes{i,3})==1
        temp2=temp;
    else
        temp2=mean(temp);
    end
    Disease_sig2(i,:)=temp2;
end

Disease_sig=Disease_sig2;
for i=1:12328
    for j=1:ns
        if Disease_sig(i,j)>1
            Disease_sig(i,j)=1;
        elseif Disease_sig(i,j)<-1
            Disease_sig(i,j)=-1;
        elseif Disease_sig(i,j)>-1&&Disease_sig(i,j)<1
            Disease_sig(i,j)=0;
        end
    end
end

%% entrez signature

load DrugSig %
nc=length(drugname);
Scores=[];
% tic
for i=1:nc
    temp=[];
        for j=1:ns
            temp(j,1)=DrugSigScore(Disease_sig(:,j),drugsigup{i},drugsigdown{i});
        end
        Scores(i,:)=temp';
end
[~,score_y]=size(Scores);
% toc
%% analysis the drugs in the results
predicted=[];predicted_value=[];
score_y1=1;
% score_y=114;
%1~75; 76~96; 97~114;115 116, patients and stages
for i=score_y1:score_y
    temp=Scores(:,i);
    [tv,tr]=sort(temp,'ascend');
    for j=1:30
        predicted{j,i}=drugname{tr(j),1};
        predicted_value(j,i)=tv(j);
    end
end
nresult1=10;
z=0;sd=[];sd_value=[];
for i=1:nresult1
    for j=score_y1:score_y
        z=z+1;
        sd{z,1}=predicted{i,j};
        sd_value(z,1)=predicted_value(i,j);
    end
end
markd=zeros(z,1);
z1=0;
for i=1:z-1
    if markd(i)==0
        z1=z1+1;
        markd(i)=z1;
        for j=i+1:z
            if markd(j)==0&&strcmp(sd{i,1},sd{j,1})
                markd(j)=z1;
            end
        end
    end
end
t1=[];t2=[];t3=[];
for i=1:z1
    pos=find(markd==i);
    t1{i,1}=sd{pos(1),1};
    t2(i,1)=length(pos);
    t3(i,1)=sum(sd_value(pos));
end
[~,t2r]=sort(t2,'descend');
resultsTOPvotes=[];
for i=1:z1
    resultsTOPvotes{i,1}=t1{t2r(i),1};
    resultsTOPvotes{i,2}=t2(t2r(i),1);
end
