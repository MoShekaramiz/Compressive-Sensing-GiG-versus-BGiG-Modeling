function [L_iter,e_Noise_iter,TAU_ITER,X_est] = Fun_GiG_DiffPrecis_VB(Y,A)
% -------------------------------------------------------------------------
% "OBL_VB" Algorithm
% This is based on the regular update rules of variables and parameters
% usig VB. It does have a greedy-based criterion.

% This code is written by Mohammad Shekaramiz
% Utah Valley University
% Last Updated: 10/19/2022
% -------------------------------------------------------------------------
% This function performs "OBL_VB" algorithm to solve the SMV problems. 
% Note: This algorithm does not accout for signals with unknown clustered 
% patterns, but rather solves for an ordinary SBL algorithm.  

% ###########
% INPUTS:
%     Y                :  The measurement matrix/vector (for SMV problem)
%     A                :  The sensing (measurement) marix
% ###########
% OUTPUTS
%    X_est             :  The estimated solution vector
% -------------------------------------------------------------------------
Phi = A' * A; % Since we will need thisa couple of times, we compute it in advance.
[ Row_A , Col_A ] = size( A );
% -------------------------------------------------------------------------
% SOLUITION COMPONENT: INITIALIZATION OF THE SOLUTION VECTOR "X".
a_0 = 1e-6;
b_0 = 1e-3;
% TAU = gamrnd( a_0 , 1 / b_0 , 1 , 1 )+1e-2;  % Initialization.
% if TAU == 0
TAU = a_0 / b_0; % end
X_BAYES = ( 1 / sqrt( TAU ) ) .* randn( Col_A , 1 );
TAU  = TAU * ones( 1, Col_A );
TAU_ITER = TAU;
%--------------------------------------------------------------------------
% Another way of initialzing the sparse solution.
% X_BAYES = A' * Y;
% Sbar = abs( X_BAYES ./ max( X_BAYES ) );
% -------------------------------------------------------------------------
% NOISE COMPONENT: INITIALIZATION OF THE ERROR PRECISION TERM "e_Noise".
theta0 = 1e-6; % 1e-6;
theta1 = 1e-6; % 1e-6;
% e_Noise = gamrnd( theta0 , 1 / theta1 , 1 )
% if e_Noise == 0
e_Noise = theta0 / theta1; %end

%--------------------------------------------------------------------------
% MAIN BODY OF THE CODE. 
% RECONSTRUCTION OF THE SPARSE SOLUTION FROM THE OBSERVED NOISY MEASUREMENTS.
%--------------------------------------------------------------------------
%--------------------------------------------------------
% Initial measure of negative-marginalized log-likelihood
Sig0_inv = 1 / e_Noise * eye( Row_A ) + A * diag(1./TAU ) *  A';% + 1e-10*eye(Row_A);
Sig0 = pinv( Sig0_inv );
Lt_0 = log( det( Sig0_inv ) ) + Y' * Sig0 * Y;
Lcheck = zeros( 1 , 2 );
Lcheck( 2 ) = Lt_0;
Like_err = 1e-8;
L_iter = [];
%--------------------------
MAX_ITER = 3000;
ITER = 0;
%----------
% Tau_iter = TAU;
e_Noise_iter = e_Noise;
%-----------------------
while (abs(Lcheck(2)-Lcheck(1))/abs(Lcheck(1))>Like_err  || norm(X_BAYES)==0)
    ITER = ITER + 1;
    if ITER == MAX_ITER
         % disp(ITER);
        break
    end
    
    Lcheck( 1 ) = Lcheck( 2 );
    L_iter = [ L_iter , Lcheck( 2 ) ];    
    %----------------------------------------------------
    % Updating the hidden variables; "X_BAYES".
    %----------------------------------------------------
    %----------------------------------
    % UPDATING THE SPARSE SOLUTION "X".
    %----------------------------------
    try
        SIGMA_BAYES = pinv( diag(TAU) + e_Noise * Phi );
    catch
    end
    X_BAYES = e_Noise * SIGMA_BAYES * A' * Y;
    if 1
        Threshold = 1e-6; % We zero out the elements of the solution with amplitude less that 1e-6;
        X_BAYES ( abs( X_BAYES ) < Threshold ) = 0;
    end
    
    % ---------------------------------------------------------------------
    % Adding the Greedy criterion.
    GVB = 0;
    if GVB
        if ITER > 20
            [ ~ , Indx_X ] = sort( abs(X_BAYES) , 'Ascend' );
            X_BAYES( Indx_X( 1 : Col_A-round(Row_A/2) ) ) = 0;
        end
    end
    % ---------------------------------------------------------------------
    
    % ---------------------------------------------------------------------
    % Updating the unknown parameters
    %-----------------------------------
    % UPDATING "TAU"; the precisions on the elements of the solution vector "X".
    %-----------------
    for kk = 1 : Col_A
        TAU(kk) = (a_0+1/2)/(b_0+0.5*((X_BAYES(kk))^2+SIGMA_BAYES(kk,kk)));
%         if TAU(kk)>100
%             TAU(kk)=1e10;
%         end
        if GVB
            if ~sign( abs( X_BAYES( kk ) ) )
                TAU( kk ) = 1e10;
            end
        end
    end
    TAU_ITER = [ TAU_ITER ; TAU ];
    %--------------------
    % NOISE COMPONENT.
    %--------------------
    Denom2 = trace( ( X_BAYES * X_BAYES' + SIGMA_BAYES ) * Phi );
    % Denom2=ones(1,Col_A)*((X_BAYES*X_BAYES'+SIGMA_BAYES).*Phi)*ones(Col_A,1);
    Denom1 = Y' * Y - 2 * X_BAYES' *A' * Y + Denom2;
    e_Noise = ( theta0 + 0.5 * Row_A ) / ( theta1 + 0.5 * Denom1 );
    e_Noise_iter = [ e_Noise_iter , e_Noise ];
    % --------------------------------------------------------------
    Sig0_inv =  1 / e_Noise * eye( Row_A ) + ( A * diag( 1./TAU ) * A' );%+1e-10*eye(Row_A);
    Sig0 = pinv( Sig0_inv );
    Lcheck( 2 ) = log( det( Sig0_inv ) ) + Y' * Sig0 * Y;
end
 
L_iter = [ L_iter , Lcheck( 2 ) ];
X_est = X_BAYES;