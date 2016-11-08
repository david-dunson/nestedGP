%Purpose:
%MCMC for nGP with m=2, n=1 and exact or approximate transitional density.

%Input:
%J: scaler, number of observations
%Y: 1*J, first value for initial value
%deltas: 1*J, duration between observations,
%             including the duration between time zero and first value
%sigma2_eps: scaler, variance of observation
%sigma2_U: scaler, variance of velocity
%sigma2_A: scaler, variance of acceleration
%nIter,BurnIn,Thin: scalers, parameters for MCMC
%exact: true/false, exact or approximate transitional density 

%Output:
%sigma2_out:(nIter-BurnIn)/Thin*3
%U_out: (nIter-BurnIn)/Thin*J
%V_out: (nIter-BurnIn)/Thin*J
%A_out: (nIter-BurnIn)/Thin*J

function [sigma2_out,U_out,V_out,A_out] = ...
    MCMC_nGP(J,Y,deltas,sigma2_eps,sigma2_U,sigma2_A,...
    nIter,BurnIn,Thin,exact)

%specify MCMC output matrixs:
sigma2_out=zeros((nIter-BurnIn)/Thin,3);
U_out=zeros((nIter-BurnIn)/Thin,J);
V_out=zeros((nIter-BurnIn)/Thin,J);
A_out=zeros((nIter-BurnIn)/Thin,J);


%state space model specifications using approximate transitional density.
%Z matrix:
m=zeros(1,3);
m(1,1)=1;
Z=ssmat(m);

%T matrix
m=eye(3);
m(1,2)=-9;
m(2,3)=-9;
dmmask=(m==-9);
dvec=repmat(deltas,2,1);
T_app=ssmat(m,[],dmmask,dvec,[]);

%R matrix
m=zeros(3,2);
m(2,1)=1;
m(3,2)=1;
R_app=ssmat(m);

%H matrix
m=sigma2_eps;
H=ssmat(m);

%Q matrix
m=zeros(2,2);
m(1,1)=-9;
m(2,2)=-9;
dmmask=(m==-9);
dvec=[sigma2_U.*deltas;sigma2_A.*deltas];
Q_app=ssmat(m,[],dmmask,dvec,[]);

nGPapp =  ssmodel(struct('type', 'nGPapp'), H, Z, T_app, R_app, Q_app);

%state space model specifications using exact transitional density.
if exact==true
    %T matrix
    m=eye(3);
    m(1,2)=-9;
    m(1,3)=-9;
    m(2,3)=-9;
    dmmask=(m==-9);
    dvec=[deltas;deltas.^2./2;deltas];
    T_ex=ssmat(m,[],dmmask,dvec,[]);
    
    %Q matrix
    m=repmat(-9,3,3);
    dmmask=(m==-9);
    dvec=[deltas.^3./3*sigma2_U+deltas.^5./20*sigma2_A;...
          deltas.^2./2*sigma2_U+deltas.^4./8*sigma2_A;...
          deltas.^3./6*sigma2_A;...
          deltas.^2./2*sigma2_U+deltas.^4./8*sigma2_A;...
          deltas*sigma2_U+deltas.^3./3*sigma2_A;...
          deltas.^2./2*sigma2_A;...
          deltas.^3./6*sigma2_A;...
          deltas.^2./2*sigma2_A;...
          deltas*sigma2_A ...
         ];
    Q_ex=ssmat(m,[],dmmask,dvec,[]);
    
    %R matrix
    m=eye(3);
    R_ex=ssmat(m);
    
    nGPex=ssmodel(struct('type', 'nGPex'), H, Z, T_ex, R_ex, Q_ex);
end

%MCMC updates:
for g=1:nIter
    if exact==false
        
        %update the states:
        [alphatilde epstilde etatilde] = simsmo(Y, nGPapp);
         
        %update the variances
       sigma2_eps = 1/gamrnd((0.01+J)/2,2/(0.01+sum(epstilde.^2)));
       sigma2_U = 1/gamrnd((0.01+J)/2,2/(0.01+sum(etatilde(1,:).^2./deltas)));
       sigma2_A = 1/gamrnd((0.01+J)/2,2/(0.01+sum(etatilde(2,:).^2./deltas)));

%         sigma2_eps = 1/gamrnd((0.1+J)/2,2/(0.1+sum(epstilde.^2)));
%         sigma2_U = 1/gamrnd((4+J)/2,2/(4+sum(etatilde(1,:).^2./deltas)));
%         sigma2_A = 1/gamrnd((4+J)/2,2/(10^12+sum(etatilde(2,:).^2./deltas)));


        
%          sigma2_eps = 1/gamrnd((2+J)/2,2/(1+sum(epstilde.^2)));
%          sigma2_U = 1/gamrnd((2+J)/2,2/(10^8+sum(etatilde(1,:).^2./deltas)));
%          sigma2_A = 1/gamrnd((2+J)/2,2/(10^8+sum(etatilde(2,:).^2./deltas)));
% 
%           sigma2_eps = 1/gamrnd((J-1)/2,2/(sum(epstilde.^2)));
%           sigma2_U = 1/gamrnd((J-1)/2,2/(sum(etatilde(1,:).^2./deltas)));
%           sigma2_A = 1/gamrnd((J-1)/2,2/(sum(etatilde(2,:).^2./deltas)));
%            
        %update SSM matrix
        %H matrix
        H=setmat(H,sigma2_eps);
        nGPapp=myset(nGPapp,'H',H);
        
        %Q matrix
        dvec=[sigma2_U.*deltas;sigma2_A.*deltas];
        Q_app=mysetdvec(Q_app,dvec);
        nGPapp=myset(nGPapp,'Q',Q_app);
        
    else %for the exact case
        %update the states:
        [alphatilde epstilde etatilde] = simsmo(Y, nGPex);
        
        %update the variances
        sigma2_eps = 1/gamrnd((0.1+J)/2,2/(0.1+sum(epstilde.^2)));
%        sigma2_eps = 1/gamrnd((J-1)/2,2/(sum(epstilde.^2)));
          
        %H matrix for nGPapp and nGPex
        H=setmat(H,sigma2_eps);
        nGPapp=myset(nGPapp,'H',H);
        nGPex=myset(nGPex,'H',H);
        
        %proposal sigma2_U_prop and sigma2_A_prop;
        [~, ~, etatilde_app] = simsmo(Y, nGPapp);
         sigma2_U_prop = 1/gamrnd((0.1+J)/2,2/(0.1+sum(etatilde_app(1,:).^2./deltas)));
         sigma2_A_prop = 1/gamrnd((0.1+J)/2,2/(0.1+sum(etatilde_app(2,:).^2./deltas)));
%          sigma2_U_prop = 1/gamrnd((J-1)/2,2/(sum(etatilde_app(1,:).^2./deltas)));
%          sigma2_A_prop = 1/gamrnd((J-1)/2,2/(sum(etatilde_app(2,:).^2./deltas)));  
          
        alpha=logpdf_omega_ex(etatilde,sigma2_U_prop,sigma2_A_prop,deltas)-...
              logpdf_omega_ex(etatilde,sigma2_U,sigma2_A,deltas)+...
              logpdf_omega_app(etatilde_app,sigma2_U,sigma2_A,deltas)-...
              logpdf_omega_app(etatilde_app,sigma2_U_prop,sigma2_A_prop,deltas);
        u=rand(1);
        if log(u)<min(0,alpha)
            sigma2_U=sigma2_U_prop;
            sigma2_A=sigma2_A_prop;
            
            %Q matrix for nGPapp
            dvec=[sigma2_U.*deltas;sigma2_A.*deltas];
            Q_app=mysetdvec(Q_app,dvec);
            nGPapp=myset(nGPapp,'Q',Q_app);
            
            %Q matrix for nGPex
            dvec=[deltas.^3./3*sigma2_U+deltas.^5./20*sigma2_A;...
                  deltas.^2./2*sigma2_U+deltas.^4./8*sigma2_A;...
                  deltas.^3./6*sigma2_A;...
                  deltas.^2./2*sigma2_U+deltas.^4./8*sigma2_A;...
                  deltas*sigma2_U+deltas.^3./3*sigma2_A;...
                  deltas.^2./2*sigma2_A;...
                  deltas.^3./6*sigma2_A;...
                  deltas.^2./2*sigma2_A;...
                  deltas*sigma2_A...
                 ];
            Q_ex=mysetdvec(Q_ex,dvec);
            nGPex=myset(nGPex,'Q',Q_ex);
        end
    end
    %save results
    if (g>BurnIn & rem(g/Thin,1)==0)
        %disp(g);
        sigma2_out((g-BurnIn)/Thin,:)=[sigma2_eps,sigma2_U,sigma2_A];
        U_out((g-BurnIn)/Thin,:,:)=alphatilde(1,:);
        V_out((g-BurnIn)/Thin,:,:)=alphatilde(2,:);
        A_out((g-BurnIn)/Thin,:)=alphatilde(3,:);
        
    end
end

