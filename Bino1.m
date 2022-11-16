function B = Bino1( N , SigmaS , gamma )


% B = nchoosek( N , SigmaS ) *  ( gamma ^ SigmaS ) * ( 1 - gamma ) ^ ( N - SigmaS );

B1 = N : -1 : N-SigmaS+1;
B1 = prod( B1 );
Sfact = factorial( SigmaS );
B = ( B1 / Sfact) * exp( SigmaS * log( gamma ) + ( N - SigmaS ) * log( 1 - gamma ) );