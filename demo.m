%plot 
    file_name=strcat('true_1','.mat');
    load(file_name);
    
    %plot(Mass/1000,Y/1000,'.b','MarkerSize',2);
    %hold on
    plot(Mass/1000,Y/1000,'-b','LineWidth',1);
    %hold off
    xlabel('km/z');
    ylabel('intensities \times 10^3');
    
%Fit 
    file_name=strcat('true_1','.mat');
    load(file_name);
    
    idx=Mass>20000 & Mass<25000;
    
    Y=Y(idx)/1000;
    deltas=[mean(diff(Mass(idx)));diff(Mass(idx))]/1000;
    J=length(Mass(idx));
    
     nIter=1500; %Trial2
     BurnIn=500;
     Thin=1;
     
     sigma2_eps=0.05;
     sigma2_U=1.4*10^6;
     sigma2_A=5;
    
    
     [sigma2_out,U_out,V_out,A_out]=...
         MCMC_nGP(J,Y',deltas',sigma2_eps,sigma2_U,sigma2_A,nIter,BurnIn,Thin,false);

 %summary  
    U_means=squeeze(mean(U_out,1));
    
    plot(Mass(idx)/1000,Y,'-b','LineWidth',1);
    xlabel('km/z');
    ylabel('intensities \times 10^3');
    hold on;
    plot(Mass(idx)/1000,U_means,'-b','LineWidth',1, 'Color', 'red');
    hold off; 
    
    
    
    


