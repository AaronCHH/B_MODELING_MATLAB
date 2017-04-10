%% 
clc;clear;

%% 
load Workshop_4_1.mat;
Z1 = Z;
load Workshop_4_2.mat;
Z2 = Z;

%%program to determine the orders of an ARMA model >>
 
pvec=[0,1,2]; % orders of AR 
qvec=[0,1,2]; % orders of MA

%% pvec and qvec are the orders of AR and MA,
%% respectively, which are going to be tested trough %% the program 


%% parameters and variables which are used in the 
%% program are presented here 

np=length(pvec);
nq=length(qvec);
aicsave=-99*ones(np,nq);
fpesave=-99*ones(np,nq);
minaic=1e+6;

%% loop to test different models 
%% using Akaike Information Criteria (AIC)

for pp=1:np 
    p=pvec(pp);
    for qq=1:nq
        q=qvec(qq);
        if p+q ~=0
           orders=[p q];
            m=armax(Z,orders);      

% m is a structure, which contains the parameters
% associated with a specific order

            aicsave(pp,qq)=aic(m);

%% aicsave, saves the AIC associated with the 
%% model, m. 

           fpesave(pp,qq)=fpe(m);
            if aicsave(pp,qq) < minaic
                minaic=aicsave(pp,qq); % save the min
                pbest=p;
                qbest=q;
                mbest=m;

%% finally, mbest saves the structure of the 
%% model with minimum AIC among the others 

            end
        end
    end
end
