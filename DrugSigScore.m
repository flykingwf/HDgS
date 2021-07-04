function [ score1 ] = DrugSigScore( diseasesig,drugup,drugdown )
%DRUGSIGSCORE Summary of this function goes here
%   Detailed explanation goes here
[x,~]=size(diseasesig);
diseaseup=find(diseasesig==1);
diseasedown=find(diseasesig==-1);
z1=length(diseaseup);z2=length(diseasedown);
z3=length(drugup);z4=length(drugdown);
score2=length(intersect(diseaseup,drugdown))+length(intersect(diseasedown,drugup));
score1=-score2/(z1+z2);

end
