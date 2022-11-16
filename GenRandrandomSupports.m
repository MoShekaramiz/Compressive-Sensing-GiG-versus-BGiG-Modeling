function S_Orig_Rand = GenRandrandomSupports( N , K )
% N: The lenght of the Signal
% K: Sparsity Level (The number of Non-zeros)
Indx = randperm(N);
S_Orig_Rand = Indx(1 : K);
S_Orig_Rand = sort( S_Orig_Rand );
% save S_Orig_Rand S_Orig_Rand