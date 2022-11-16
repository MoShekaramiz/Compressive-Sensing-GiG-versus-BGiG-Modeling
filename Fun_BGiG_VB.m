function [ L_iter, Epsilon_iter , Tau_iter , Gama_iter , Sbar_Collect , S_est , X_est ] = Fun_BGiG_VB( Y, A )
% -------------------------------------------------------------------------
% "O-SBL_VB" Algorithm.
% It has some issues with the selection of hyper-parameters on gamma
% distribtion for small sampling rates.
% This algorithm is based on the regular update rules of variables and 
% parameters usig VB. 
% It does not have any greedy-based criterion as we have in G-OSBL(VB).
% Note: This algorithm does not accout for signals with unknown clustered 
% patterns, and only solves for an ordinary SBL algorithm.  

% This code is written by Mohammad Shekaramiz
% Utah Valley University
% Last Updated: 10/03/2022
% -------------------------------------------------------------------------
% This function performs "OSBL_VB" algorithm to solve the SMV problem. 

% ###########
% INPUTS:
%     Y                :  The measurement vector (for SMV problem)
%     A                :  The sensing (measurement) marix
% ###########
% OUTPUTS
%    X_est             :  The estimated solution vector
%    S_est             :  The estimated support vector
%    Sbar_Collect      :  The collection of Support learning vector vs. iterations.
%    L_iter            :  The collection of negative marginalized-log-likelihood vs. iterations.
%    Epsilon_iter      :  The collection of the estimated precision on the measurement noise vs. iterations.
%    Tau_iter          :  The collection of the estimated precision on the components of the solution vs. iterations.
%    Gama_iter         :  The collection of the coefficient parameters inluencing on the Support learning vector vs. iterations.
% -------------------------------------------------------------------------
Phi = A' * A; % Since we will need this a couple of times, we compute it in advance.
[ Row_A , Col_A ] = size( A );
% -------------------------------------------------------------------------
% Initialization of the hidden variables and parameters of the model.
% -------------------------------------------------------------------------
% SPARSE SUPPORT COMPONENT: INITIALIZATION OF VECTOR "S".
Alpha0_S = 0.01; Beta0_S = 0.99; % Case 1
% Alpha0_S = 0.1; Beta0_S = 0.9; % Case 2
% Alpha0_S = 1.4; Beta0_S = 2.00; % Case 3
Gama = betarnd( Alpha0_S , Beta0_S , Col_A , 1 ); 
S = Bernoulli( Gama , 1 );
while ~any( S )
    Gama = betarnd( Alpha0_S , Beta0_S , Col_A , 1 );
    S = Bernoulli( Gama , 1 );
end
Sbar = double( S );
% -------------------------------------------------------------------------
% SOLUITION COMPONENT: INITIALIZATION OF THE SOLUTION VECTOR "X".
a_0 = 1e-3;
b_0 = 1e-3;
% TAU = gamrnd( a_0 , 1 / b_0 , 1 , 1 )+1e-3;  % Initialization.
TAU = a_0 / b_0; % end
X_BAYES = ( 1 / sqrt( TAU ) ) .* randn( Col_A , 1 );
X_BAYES = X_BAYES .* Sbar; % Generating an initial random sparse solution.
%--------------------------------------------------------------------------
% Another way of generating the initial guess on the solution.
if 1
    X_BAYES = A' * Y;
    Sbar = abs( X_BAYES ./ max( abs( X_BAYES ) ) );
end
% -------------------------------------------------------------------------
% NOISE COMPONENT: INITIALIZATION OF THE ERROR PRECISION TERM "Epsilon".
theta0 = 1e-6;
theta1 = 1e-6;
% Epsilon = gamrnd( theta0 , 1 / theta1 , 1 )
Epsilon = theta0 / theta1;
% -------------------------------------------------------------------------
% Initial measure of negative-marginalized log-likelihood
Sig0_inv = 1/Epsilon*eye(Row_A) + 1/TAU*(A*diag(Sbar.^2)*A');% + 1e-10*eye(Row_A);
Sig0 = pinv( Sig0_inv );
Lt_0 = log( det( Sig0_inv ) ) + Y' * Sig0 * Y;
Lcheck = zeros( 1 , 2 );
Lcheck( 2 ) = Lt_0;
Like_err = 1e-8;
L_iter = [];
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% MAIN BODY OF THE CODE. 
% RECONSTRUCTION OF THE SPARSE SOLUTION FROM THE OBSERVED NOISY MEASUREMENTS.
% -------------------------------------------------------------------------
MAX_ITER = 500;
ITER = 0;
% --------------------------------
Sbar_Collect = Sbar;
Tau_iter = TAU;
Epsilon_iter = Epsilon;
Gama_iter = Gama;
Alpha1_S = Alpha0_S + Sbar;
Beta1_S = Beta0_S + ( 1 - Sbar ); 
Threshold = 1e-6;
%----------------------------------

while (abs(Lcheck(2)-Lcheck(1))/abs(Lcheck(1))>Like_err  || norm(X_BAYES)==0) && ( ITER < MAX_ITER )
    ITER = ITER + 1;
    Lcheck( 1 ) = Lcheck( 2 );
    L_iter = [ L_iter , Lcheck( 2 ) ];
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % Updating the hidden variables; "X_BAYES and Sbar".
    %----------------------------------------------------
    Sbar_comptStuff = Sbar * Sbar' + diag( Sbar .*( 1 - Sbar ) );
    %----------------------------------
    % UPDATING THE SPARSE SOLUTION "X".
    %----------------------------------
    SIGMA_BAYES = pinv( TAU*eye(Col_A)+Epsilon*( Phi .*Sbar_comptStuff ) );
    X_BAYES = Epsilon * SIGMA_BAYES * diag( Sbar ) * A' * Y;
    if 1
        % We zero out the elements of the solution with amplitude less that 1e-6;
        X_BAYES ( abs( X_BAYES ) < Threshold ) = 0;
    end
    % ---------------------------------------------------------------------
    % UPDATING THE ELEMENTS OF THE SUPPORT VECTOR "S".
    % ---------------------------------------------------------------------
    for n = 1 : Col_A
        %------------------  
        % UPDATING Y_k1had.
        %------------------
        Y_K_HAD = Y - A *(Sbar.*X_BAYES) + A(:,n) *(Sbar(n).* X_BAYES(n));
        %-------------------
        C_n = exp( psi( Beta1_S( n ) ) - psi( Alpha1_S( n ) ) );
        Kappa_n = exp(0.5*Epsilon*((X_BAYES(n)^2+SIGMA_BAYES(n,n))*(norm(A(:,n))^2)-2*X_BAYES(n)*(Y_K_HAD'*A(:,n) )));
        Q_denom = 1 / ( 1 + C_n * Kappa_n );
        Q_denom = min( Q_denom , 1e100 );  % To avoid numerical issues.
        Q_denom = max( Q_denom , 1e-100 ); % To avoid numerical issues.
        Sbar ( n ) = Q_denom;   
    end
    Sbar_Collect = [ Sbar_Collect , Sbar ];
    % ---------------------------------------------------------------------
    % ---------------------------------------------------------------------
    % Updating the unknown parameters
    % ---------------------------------------------------------------------
    % UPDATING "TAU"; the precision on the solution vector "X".
    % ---------------------------------------------------------------------
    TAU = (a_0+0.5*Col_A)/(b_0+0.5*(norm(X_BAYES,2)^2+trace(SIGMA_BAYES)));
    Tau_iter = [ Tau_iter , TAU ];
    %----------------------------------------------------------------------
    % UPDATING "Epsilon"; the precision on the measurement noise.
    % ---------------------------------------------------------------------
    Sbar_comptStuff = Sbar * Sbar' + diag( Sbar .*( 1 - Sbar ) );
    % Denom2=trace((X_BAYES*X_BAYES'+SIGMA_BAYES) * (Phi.* Sbar_comptStuff));
    Denom2=ones(1,Col_A)*((X_BAYES*X_BAYES'+SIGMA_BAYES).*(Phi).*Sbar_comptStuff)*ones(Col_A,1); % This is faster.
    Denom1 = Y' * Y - 2 * ( X_BAYES .* Sbar )' *A' * Y + Denom2;
    Epsilon = ( theta0 + 0.5 * Row_A ) / ( theta1 + 0.5 * Denom1 );
    Epsilon_iter = [ Epsilon_iter , Epsilon ];
    % ---------------------------------------------------------------------
    % UPDATING "Gama" WHICH IS USED FOR MAKING DESICION ON THE 
    % ENTRIES OF THE SUPPORT VECTOR "S".
    % --------------------------------------------------------------------- 
    Alpha1_S = Alpha0_S + Sbar;
    Beta1_S = Beta0_S + ( 1 - Sbar );
    % We compute Gamma only for investigating the learning process. It is
    % not required to be computed in the algorithm.
    Gama = Alpha1_S ./ ( Alpha1_S + Beta1_S );
    Gama_iter = [ Gama_iter , Gama ];
    % --------------------------------------------------------------------- 
    % --------------------------------------------------------------------- 
    % Computing the Negative log-likelihood for the stopping condition.
    Sig0_inv = 1/Epsilon*eye(Row_A)+1/TAU*(A*diag(Sbar.^2)*A');%+1e-10*eye(Row_A);
    Sig0 = pinv( Sig0_inv );
    Lcheck( 2 ) = log( det( Sig0_inv ) ) + Y' * Sig0 * Y;
end
 
L_iter = [ L_iter , Lcheck( 2 ) ];
S_est = Sbar;
X_est = X_BAYES;
