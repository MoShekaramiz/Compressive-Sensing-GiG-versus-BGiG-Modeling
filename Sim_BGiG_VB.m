clc, clear, close all
SNR = [ 5 , 10 , 15 , 20 , 25 ];
noise_var = 1 ./ (10 .^( SNR / 10 ) );
noise_var = noise_var( 5 );
SNR = SNR( 5 );
disp( SNR );
%------------------------------------------------
% GENERATING MEASUREMENTS AND THE TRUE SOLUTION.
% M : Number of measurements
Case = 1;
M = 80;
N = 100;
Ksparsity = 25;
Threshol = 1e-3;
Eval_Tot = 0;
Max_iter = 1; % 1000
for ITER = 1 : Max_iter
    %disp(ITER);
     
    A = randn( M , N ) / sqrt( M );
    A = normc( A );
    % S_Orig = GenRandClusteredSupports( N , Ksparsity1 );
    S_Orig = GenRandrandomSupports( N , Ksparsity );
    X_TRUE = zeros( N , 1 );
    X_TRUE( S_Orig ) = randn( Ksparsity , 1 );
    S_True = zeros( size( X_TRUE ) );
    S_True (S_Orig ) = 1;
    NOISE = sqrt(noise_var) * randn( M , 1 );
    Y = A * X_TRUE + NOISE;
    %%%%%%%%%%%
    
    %-------------------------------------
    [ L_iter , e_Noise_iter , Tau_iter , Gama_iter , Sbar_Collect , S_OSBL_VB , X_OSBL_VB ] = Fun_BGiG_VB( Y , A );
    SUPP_OSBL_VB = find( S_OSBL_VB >= Threshol );
    X_OSBL_VB( abs( X_OSBL_VB ) < Threshol ) = 0 ;
    X_ERR_OSBL_VB = ERR_RATE( X_TRUE , X_OSBL_VB );
    S_OSBL_VB( S_OSBL_VB <= Threshol ) = 0;
    [ PD_OSBL_VB , PFA_OSBL_VB ] = DET_FALSE_RATE( Ksparsity, N, S_Orig, SUPP_OSBL_VB );
    Eval_OSBL_VB = [ PD_OSBL_VB, PFA_OSBL_VB,PD_OSBL_VB-PFA_OSBL_VB, X_ERR_OSBL_VB, M/N ];     
    % display( Eval_OSBL_VB );
    Eval_Tot = Eval_Tot + Eval_OSBL_VB;
end
Eval_Tot = Eval_Tot / ITER;
%save(['Eval_Tot_',num2str(M)], 'Eval_Tot');
disp( Eval_Tot );
% save(['L_',num2str(M),'_SNR_',num2str(SNR)], 'L_iter');
% save(['e__',num2str(M),'_SNR_',num2str(SNR)], 'e_Noise_iter');
% save(['Tau_',num2str(M),'_SNR_',num2str(SNR)], 'Tau_iter');
%-------------------------------------
% PLOTS
if 1
    figure; plot(X_TRUE,'b*-.','LineWidth',1.4);grid; hold on; plot(X_OSBL_VB,'r','LineWidth',1.4); 
    set(gca,'FontSize',16);
    h1 = legend('$\mathbf{x}$','$\hat{\mathbf{x}}$','Location','best'); set(h1,'FontSize',14,'Interpreter','latex');
%     filename = ['X_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.fig']; saveas(gcf,filename);
%     filename = ['X_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.bmp']; saveas(gcf,filename);



    figure; plot(Y,'b-.*','LineWidth',1.4);grid; hold on; plot(A*(X_OSBL_VB),'r','LineWidth',1.4); 
    set(gca,'FontSize',16);
    h2 = legend('$\mathbf{y}$','$\hat{\mathbf{y}}$','Location','best'); set(h2,'FontSize',14,'Interpreter','latex');
%     filename = ['Y_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.fig']; saveas(gcf,filename);
%     filename = ['Y_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.bmp']; saveas(gcf,filename);
    
    Len_Iter = length(Tau_iter);
    figure; plot(0:Len_Iter-1,ones(1,length(Tau_iter)),'b','LineWidth',2); grid; hold on;
    plot(0:Len_Iter-1,Tau_iter,'r','LineWidth',2); title('Estimate of Solution Precision', 'FontSize', 14);
    xlabel('Iteration', 'FontSize', 14); set(gca,'FontSize',16);
    h3 = legend('$\tau$','$\hat{\tau}$','Location','best'); set(h3,'FontSize', 14,'Interpreter','latex');
%     filename = ['Tau_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.fig']; saveas(gcf,filename);
%     filename = ['Tau_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.bmp']; saveas(gcf,filename);
    
    figure; plot(0:Len_Iter-1,1/noise_var*ones(1,length(Tau_iter)),'b','LineWidth',2); grid; hold on;
    plot(0:Len_Iter-1,e_Noise_iter,'r','LineWidth',2); title('Estimate of Noise Precision', 'FontSize', 14);
    xlabel('Iteration', 'FontSize', 14);  set(gca,'FontSize',16);
    h6 = legend('$\varepsilon$','$\hat{\varepsilon}$','Location','best'); set(h6,'FontSize',14,'Interpreter','latex');
%     filename = ['Noise_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.fig']; saveas(gcf,filename);
%     filename = ['Noise_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.bmp']; saveas(gcf,filename);

    figure; plot(0:Len_Iter-1,L_iter,'LineWidth',2); grid; title('Log-Marginalized-Likelihood', 'FontSize', 14);
    xlabel('Iteration', 'FontSize', 14);ylabel('$L$', 'FontSize', 14,'Interpreter','latex'); set(gca,'FontSize',16);
%     filename = ['L_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.fig']; saveas(gcf,filename);
%     filename = ['L_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.bmp']; saveas(gcf,filename); 

    figure; stem(S_True,'b*-.','LineWidth',1.4);grid; hold on; stem(S_OSBL_VB,'r','LineWidth',1.4); 
    set(gca,'FontSize',16);
    h4 = legend('$\mathbf{s}$','$\hat{\mathbf{s}}$','Location','best'); set(h4,'FontSize',14,'Interpreter','latex');
%     filename = ['S_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.fig']; saveas(gcf,filename);
%     filename = ['S_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.bmp']; saveas(gcf,filename);

    figure; imagesc(Sbar_Collect); set(gca,'FontSize',16); 
    xlabel('Iteration', 'FontSize', 14); ylabel('$\hat{\mathbf{s}}_{collect}$', 'FontSize', 24,'Interpreter','latex');
%     filename = ['Scoll_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.fig']; saveas(gcf,filename);
%     filename = ['Scoll_SNR',num2str(SNR),'_Case',num2str(Case),'_',num2str(M),'.bmp']; saveas(gcf,filename);
end




