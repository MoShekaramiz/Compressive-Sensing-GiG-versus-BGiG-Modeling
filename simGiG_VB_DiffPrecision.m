clc, clear, close all
SNR = [ 5 , 10 , 15 , 20 , 25 ];
n_Var = 1 ./ (10 .^( SNR / 10 ) );
VAR = n_Var( end );
%------------------------------------------------
% M : Number of measurements
M = 80;
N = 100;
Threshol = 1e-3;
Threshol1 = 5e-2;
Num_sim  = 1; %1000
Eval = 0; 
    for MM = 1 : Num_sim
        %disp( MM );
        Ksparsity = 25;
        
            A = randn( M , N ) / sqrt( M );
            A = normc( A );
            S_Orig = GenRandrandomSupports( N , Ksparsity );
            Ksparsity = length( S_Orig );
            X_TRUE = zeros( N , 1 );
            X_TRUE( S_Orig ) = randn( Ksparsity , 1 );
            X_TRUE(abs(X_TRUE) < Threshol1 & abs(X_TRUE) > 0) = Threshol1;
            
            NOISE = sqrt(VAR) * randn( M , 1 );
            Y = A * X_TRUE + NOISE;

        S_True = zeros( size( X_TRUE ) );
            S_True (S_Orig ) = 1;
        [L_ITER,Noise_ITER, TAU_ITER, X_OSBL_VB] = Fun_GiG_DiffPrecis_VB( Y , A );
        SUPP_OSBL_VB = find( abs(X_OSBL_VB) > Threshol );
        X_OSBL_VB( abs( X_OSBL_VB ) <= Threshol ) = 0 ;
        X_ERR_OSBL_VB = ERR_RATE( X_TRUE , X_OSBL_VB );
        [ PD_OSBL_VB , PFA_OSBL_VB ] = DET_FALSE_RATE( Ksparsity, N, S_Orig, SUPP_OSBL_VB );
        PdPfa_Diff = PD_OSBL_VB - PFA_OSBL_VB;
        Eval_OSBL_VB = [ PD_OSBL_VB, PFA_OSBL_VB, PdPfa_Diff, X_ERR_OSBL_VB , M/N ];
        % disp(Eval_OSBL_VB);
        Eval = Eval + Eval_OSBL_VB;
    end
    Eval_OBL = Eval / MM;
    disp( Eval_OBL );
    plot(X_TRUE); hold on; grid; plot(X_OSBL_VB); legend('X_{TRUE}','X_{EST}')
    figure; plot(ones(1,N)); hold on; grid; plot(TAU_ITER); legend('Tau','Tau_{est}')%(end,:)
    figure; plot(L_ITER); grid; legend('Marginalied Log-Likelihood')
    figure; plot(1/VAR*ones(1,N)); hold on; grid; plot(Noise_ITER); legend('Noise Precision','Estimated Noise Precision')
   % save(['Eval_OBL_',num2str(M)],'Eval_OBL');