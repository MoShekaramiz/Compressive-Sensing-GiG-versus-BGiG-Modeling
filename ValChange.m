function    SigmaDelta = ValChange( S )

% LocChanges = find(diff(S)) + 1;
% SigmaDelta = length( LocChanges );
SigmaDelta = sum( abs( diff( S ) ) );