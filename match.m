

% nm and tm in the seperate files
% marks=ones(20291,1);
% for i=1:12328
%     temp=[];
%     z=0;
%     maps{i,1}=LINCSgeneinfo{i,1};
%     for j=1:20291
%         if marks(j)&&LINCSgeneinfo{i,1}==probeToLINCSentrez{j,2}
%             z=z+1;
%             temp{z}=probeToLINCSentrez{j,1};
%         end
%     end
%     maps{i,2}=temp;
% end
% 
% for i=1:12328
%     if ~isempty(maps{i,2})
%         temp=maps{i,2};
%         temp2=[];
%         for j=1:length(temp)
%             for jj=1:22277
%                 if strcmp(temp{j},probe96{jj,1})
%                     temp2=[temp2,jj];
%                 end
%             end
%         end
%         maps{i,3}=temp2;
%     end
% end
% 
% for i=1:2883
%     temp=CORUM_com_F0{i,4};
%     temp2=[];
%     for j=1:length(temp)
%         temp2=[temp2,maps{temp(j),3}];
%     end
%     maps_corumF0_probes{i,1}=temp2;
% end
% save maps_corumF0_probes.mat maps_corumF0_probes
%disease mat1
load maps_corumF0_probes
% for i=1:2624
%     temp=mat1(maps_corumF1_probes{i,1},:);
%     temp2=mean(temp);
%     mat2(i,:)=temp2;
% end
% Disease_sig1=mat2(:,ns+1:ns*2)-mat2(:,1:ns);
% %disease mat2
Disease_sig2=[];
mat1=nm-tm;
[~,ns]=size(mat1);
for i=1:2883
    temp=mat1(maps_corumF0_probes{i,1},:);
    if length(maps_corumF0_probes{i,1})==1
        temp2=temp;
    else
        temp2=mean(temp);
    end
    Disease_sig2(i,:)=temp2;
end
Disease_sig=[];
for i=1:ns
    temp1=Disease_sig2(:,i);
    [tv,tr]=sort(temp1,'descend');
    temp2=[];
    for j=1:2883
        temp2(tr(j),1)=j;
    end
    Disease_sig(:,i)=temp2;
end

% load('maps_LINCS_probes.mat')
% load CORUM_com_F0
% for i=1:12328
%     temp=mat1(maps_LINCS_probes{i,3},:);
%     if length(maps_LINCS_probes{i,3})==1
%         temp2=temp;
%     else
%         temp2=mean(temp);
%     end
%     Disease_sig2(i,:)=temp2;
% end
% CORUM_DIST_mat=[];
% z=0;
% for i=1:2883
%     temp2=CORUM_com_F0{i,4};
%     if length(temp2)>2
%         temp3=Disease_sig2(temp2,:);
%         ntemp=length(temp2);
%         commat=dist(temp3); %corrcoef: pearson, corr: spearman
%         CORUM_DIST_mat{i,1}=commat;
%     end
% end
% CORUM_DIST_vector=[];
% for i=1:2883
%     ctemp=CORUM_DIST_mat{i,1};
%     if length(CORUM_com_F0{i,4})==1
%         CORUM_DIST_vector{i,1}=1;
%     elseif length(CORUM_com_F0{i,4})==2
%         CORUM_DIST_vector{i,1}=[0.5,0.5];
%     else
%         cv=sum(ctemp);
%         cm=mean(cv);
%         ct2=cm./cv;
%         cv2=ct2/sum(ct2);
%         CORUM_DIST_vector{i,1}=cv2;
%     end
% end
% mat1=FromLincsToCORUM_PCCweighted(Disease_sig2,CORUM_com_F0,CORUM_DIST_vector);
% 
% CORUM_PCC_mat=[];
% for i=1:2883
%     temp2=CORUM_com_F0{i,4};
%     if length(temp2)>2
%         temp3=Disease_sig2(temp2,:);
%         ntemp=length(temp2);
%         commat=corr(temp3','Type','Pearson'); %corrcoef: pearson, corr: spearman
%         CORUM_PCC_mat{i,1}=commat;
%     end
% end
% CORUM_PCC_vector=[];z1=0;z2=0;z3=0;
% wall=[];
% for i=1:2883
%     temp2=CORUM_com_F0{i,4};
%     nss=length(temp2);
%     if nss==1
%         wt=[1];z1=z1+1;wall=[wall,wt];
%         CORUM_PCC_vector{i,1}=wt;
%     elseif nss==2
%         wt=[0.5,0.5];z2=z2+1;wall=[wall,wt];
%         CORUM_PCC_vector{i,1}=wt;
%     else
%         commat=CORUM_PCC_mat{i,1};
%         comsum=sum(commat)-1;
%         wt=comsum/sum(comsum);
%         z3=z3+1;wall=[wall,wt];
%         CORUM_PCC_vector{i,1}=wt;
%     end
% end
% mat2=FromLincsToCORUM_PCCweighted(Disease_sig2,CORUM_com_F0,CORUM_PCC_vector);
% 
% CORUM_Spear_mat=[];
% for i=1:2883
%     temp2=CORUM_com_F0{i,4};
%     if length(temp2)>2
%         temp3=Disease_sig2(temp2,:);
%         ntemp=length(temp2);
%         commat=corr(temp3','Type','Spearman'); %corrcoef: pearson, corr: spearman
%         CORUM_Spear_mat{i,1}=commat;
%     end
% end
% CORUM_Spear_vector=[];z1=0;z2=0;z3=0;
% wall=[];
% for i=1:2883
%     temp2=CORUM_com_F0{i,4};
%     nss=length(temp2);
%     if nss==1
%         wt=[1];z1=z1+1;wall=[wall,wt];
%         CORUM_Spear_vector{i,1}=wt;
%     elseif nss==2
%         wt=[0.5,0.5];z2=z2+1;wall=[wall,wt];
%         CORUM_Spear_vector{i,1}=wt;
%     else
%         commat=CORUM_Spear_mat{i,1};
%         comsum=sum(commat)-1;
%         wt=comsum/sum(comsum);
%         z3=z3+1;wall=[wall,wt];
%         CORUM_Spear_vector{i,1}=wt;
%     end
% end
% mat3=FromLincsToCORUM_PCCweighted(Disease_sig2,CORUM_com_F0,CORUM_Spear_vector);
% DS1=[];
% for i=1:ns
%     temp1=mat1(:,i);
%     [tv,tr]=sort(temp1,'descend');
%     temp2=[];
%     for j=1:2883
%         temp2(tr(j),1)=j;
%     end
%     DS1(:,i)=temp2;
% end
% DS2=[];
% for i=1:ns
%     temp1=mat2(:,i);
%     [tv,tr]=sort(temp1,'descend');
%     temp2=[];
%     for j=1:2883
%         temp2(tr(j),1)=j;
%     end
%     DS2(:,i)=temp2;
% end
% DS3=[];
% for i=1:ns
%     temp1=mat3(:,i);
%     [tv,tr]=sort(temp1,'descend');
%     temp2=[];
%     for j=1:2883
%         temp2(tr(j),1)=j;
%     end
%     DS3(:,i)=temp2;
% end

% a=sum(sum(abs(Disease_sig1-Disease_sig2)));
%% 
% load SIG_compound_DIST_F0 %DIST PCC Spear Unweighted
% load SIG_compound_PCC_F0
% load SIG_compound_Spear_F0
% load SIG_compound_Unweighted_F0
% % load SIG_compound_PCCweighted2_F1_0810
% % load SIG_compound_cellline_PCCweighted
% load FILE_level5_compound_cells
% load FILE_level5_compound_drugbank
% nc=length(SIG_compound_DIST_F0);
% SIG_compound_DIST_F0{1,1}=[];SIG_compound_DIST_F0{1,2}=[];
% 
% % Scores=zeros(nc,ns);
% Scores=[];
% for i=1:nc
%     temp=[];
%     if FILE_level5_compound_drugbank(i)%~isempty(SIG_compound{i,1})
%         for j=1:ns
%             temp(j,1)=Score_CMap(SIG_compound_DIST_F0{i,2},Disease_sig(:,j));
%         end
%         Scores(i,:)=temp';
%     end
% end
% [~,score_y]=size(Scores);
% load smallcompounds
% 
% 
% %% analysis the drugs in the results
% predicted=[];predicted_value=[];
% for i=1:score_y
%     temp=Scores(:,i);
%     [tv,tr]=sort(temp,'ascend');
%     for j=1:30
%         predicted{j,i}=smallcompounds{tr(j),1};
%         predicted_value(j,i)=tv(j);
%     end
% end
% nresult1=10;
% z=0;sd=[];sd_value=[];
% for i=1:nresult1
%     for j=1:score_y
%         z=z+1;
%         sd{z,1}=predicted{i,j};
%         sd_value(z,1)=predicted_value(i,j);
%     end
% end
% markd=zeros(z,1);
% z1=0;
% for i=1:z-1
%     if markd(i)==0
%         z1=z1+1;
%         markd(i)=z1;
%         for j=i+1:z
%             if markd(j)==0&&strcmp(sd{i,1},sd{j,1})
%                 markd(j)=z1;
%             end
%         end
%     end
% end
% t1=[];t2=[];t3=[];
% for i=1:z1
%     pos=find(markd==i);
%     t1{i,1}=sd{pos(1),1};
%     t2(i,1)=length(pos);
%     t3(i,1)=sum(sd_value(pos));
% end
% [~,t2r]=sort(t2,'descend');
% resultsTOPvotes=[];
% for i=1:z1
%     resultsTOPvotes{i,1}=t1{t2r(i),1};
%     resultsTOPvotes{i,2}=t2(t2r(i),1);
% end
% [~,t3r]=sort(t3,'ascend');
% resultsTOPvalues=[];
% for i=1:z1
%     resultsTOPvalues{i,1}=t1{t3r(i),1};
%     resultsTOPvalues{i,2}=t3(t3r(i),1);
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compound cell line signature
% load SIG_compound_DIST_F0 %DIST PCC Spear Unweighted
load SIG_compound_PCC_F0_200
% load SIG_compound_Spear_F0
% load SIG_compound_Unweighted_F0
% load SIG_compound_cellline_DIST_F0 %DIST PCC Spear Unweighted
% load SIG_compound_cellline_PCC_F0
% load SIG_compound_cellline_Spear_F0
% load SIG_compound_cellline_Unweighted_F0
% load SIG_compound_cellline
% load FILE_level5_compound_cells
load FILE_level5_compound
load FILE_level5_compound_drugbank
nc=length(SIG_compound_Unweighted_F0);
Scores=zeros(nc,ns);
% Scores=[];
s1=60;
for i=1:nc
    if FILE_level5_compound_drugbank(i)
    temp=[];
        for j=1:ns
            sigtemp=SIG_compound_Unweighted_F0{i,3};
            [x1,~]=size(sigtemp);
            if x1<s1
                sigTEMP=sigtemp;
            else
                sigTEMP=sigtemp(1:s1,:);
            end
            temp(j,1)=Score_CMap(sigTEMP,Disease_sig(:,j));
        end
        Scores(i,:)=temp';
    end
end
[~,score_y]=size(Scores);

%% analysis the drugs in the results
load knowndrugs
predicted=[];predicted_value=[];
score_y1=1;
% score_y=114;
%1~75; 76~96; 97~114;115 116
for i=score_y1:score_y
    temp=Scores(:,i);
    [tv,tr]=sort(temp,'ascend');
    for j=1:30
        predicted{j,i}=FILE_level5_compound{tr(j),1};
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
%     if i<21
%     for j=1:length(knowndrugs_lung)
%         if strcmpi(resultsTOPvotes{i,1},knowndrugs_lung{j,1})
%             resultsTOPvotes{i,1}
%         end
%     end
%     end
end

%%%%%%%%%%%%%%%%%%%%
resultsTOPvotes=[];
for i=1:z1
    resultsTOPvotes{i,1}=t1{t2r(i),1};
    resultsTOPvotes{i,2}=t2(t2r(i),1);
    if i<11
    for j=1:length(knowndrugs_colon) % breast lung colon prostate
        if strcmpi(resultsTOPvotes{i,1},knowndrugs_colon{j,1})
            resultsTOPvotes{i,:}
        end
    end
    end
end








%% analysis the drugs in the results

predicted=[];
for i=1:score_y
    temp=Scores(:,i);
    [tv,tr]=sort(temp,'descend');
    for j=1:30
        predicted{j,i}=FILE_level5_compound{tr(j),1};
    end
end

nresult1=10;
z=0;sd=[];
for i=1:nresult1
    for j=1:score_y
        z=z+1;
        sd{z,1}=predicted{i,j};
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
t1=[];t2=[];
for i=1:z1
    pos=find(markd==i);
    t1{i,1}=sd{pos(1),1};
    t2(i,1)=length(pos);
end
[~,t2r]=sort(t2,'descend');
resultsBOT=[];
for i=1:z1
    resultsBOT{i,1}=t1{t2r(i),1};
    resultsBOT{i,2}=t2(t2r(i),1);
end
