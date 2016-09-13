function [I, s, p, q, r, c, rh, ch, B, u, v] = mondriaan(A, NumProcessors, Imbalance, Permutation, Symmetry, SplitStrategy)
% mondriaan       Uses the Mondriaan algorithm to assign the nonzeros
%                 of a given sparse matrix A to NumProcessors processors
%                 such that the required communication for the matrix-
%                 vector multiplication u = A*v is as small as possible,
%         while ensuring that the number of nonzeros assigned to each
%         processor has a maximum imbalance given by Imbalance.
%
%         If desired, a permutation of A may also be specified in the
%         variable Permutation: 0 = no permutation,
%                               1 = move cut rows/columns to the front
%                                   (reverse Bordered Block Diagonal)
%                               2 = move cut rows/columns to the middle
%                                   (Separated Block Diagonal)
%                               3 = move cut rows/columns to the back.
%                                   (Bordered Block Diagonal).
%
%         If desired, a symmetry argument can be passed, affecting how 
%         the input matrix is translated from Matlab to the Mondriaan
%         library. The variable Symmetry can be set to:
%             0 = not symmetric (matrix is passed as-is),
%             1 = matrix is symmetric, (pass only the lower triangular part),
%             2 = matrix is structurally symmetric,
%                 (pass only the lower triangular part, but reconstruct
%                  correct matrix if Permutation>0)
%         
%         
%         Note that these last two input parameters are optional, while 
%         the first three input parameters are not. By default,
%         Permutation = Symmetry = 0.
%         
%         Return values:
%             I = Copy of the matrix A where the numerical value of each
%                 nonzero is replaced by the index of the processor to
%                 which this nonzero was assigned to by the Mondriaan
%                 algorithm.
%             s = The statistics vector with information about the
%                 partitioning: duration, imbalance, max communication and
%                 communication volume. If no vector distribution was
%                 requested, the latter two values will be set to -1.
%             p = Permutation vector corresponding to the matrix P (see above),
%                 so that PA = A(p,:).
%             q = Permutation vector corresponding to the matrix Q, so that
%                 AQ = A(:,q) and PAQ is given by A(p,q).
%             r = Vector containing the indices of each block visible after 
%                 permutation (thus each block in PAQ); block i starts at the
%                 row index r(i), and ends at r(i+1). Start index is inclusive,
%                 end indices are exclusive. r(1) is always 1, r(end) is always
%                 size(A,1)+1.
%             c = Vector containing the indices of each block visible after
%                 permutation; block i starts at the column index c(i) and ends
%                 at c(i+1). See above.
%             rh= The row-wise block hierarchy. For each block i, rh(i) gives
%                 its parent block. This builds a binary tree (which need not be
%                 complete), where the largest separator block (spanning all
%                 rows) is the root. Separator blocks are internal nodes of this
%                 tree; the pure nonseparated blocks are the leaf nodes.
%             ch= The column-wise block hierarchy. For each block i, rc(i) gives
%                 its column-wise parent block. See above.
%             B = A (permuted if Permutation>0) copy of the matrix A: B=PAQ.
%             v = Processor indices to which the components of the vector v
%                 are assigned. These are never permuted so to correspond with I.
%             u = Processor indices to which the components of the vector u
%                 are assigned. These are never permuted.
%         
%         Usage:
%             [I, s] = mondriaan(A, NumProcessors, Imbalance, Permutation, Symmetry);
%             [I, s, p, q] = mondriaan(A, NumProcessors, Imbalance, Permutation, Symmetry);
%             [I, s, p, q, r, c, rh, ch, B] = mondriaan(A, NumProcessors, Imbalance, Permutation, Symmetry);
%             [I, s, p, q, r, c, rh, ch, B, u, v] = mondriaan(A, NumProcessors, Imbalance, Permutation, Symmetry);
%

if (nargin < 3)
    disp 'Error: first three arguments to mondriaan function are mandatory!';
    return;
end
if (NumProcessors < 0)
    disp('Invalid value for NumProcessors');
    return;
end
if (Imbalance <= 0)
    disp('Invalid value for Imbalance');
    return;
end
if (issparse(A) ~= 1)
    disp('Input matrix A is not sparse');
    return;
end
if (nargin < 4)
    Permutation = 0;
elseif( Permutation < 0 || Permutation > 3)
    disp('Invalid value for Permutation');
    return;
end
if (nargin < 5)
    Symmetry = 0;
elseif( Symmetry < 0 || Symmetry > 2)
    disp('Invalid value for Symmetry');
    return;
elseif( Symmetry > 0 )
    disp('------------------------------------------------------------------------------------------------');
    disp('|Warning: the input matrix is assumed symmetric, this call will not fail if it actually is not!|');
    disp('------------------------------------------------------------------------------------------------');
    disp('Info: recommend symfinegrain split strategy');
    disp('Info: forcing symmetry overrides "UseSingleEntry"  option to "yes"');
    disp('Info: forcing symmetry overrides "SingleEntryType" option to "lower"');
end
if (nargin<6)
    SplitStrategy = -1;
elseif (SplitStrategy <-1 || SplitStrategy > 8)
    disp('Invalid value for SplitStrategy, falling back to the default one.')
    SplitStrategy = -1;
end

if (Symmetry > 0)
    Pass = tril(A);
else
    Pass = A;
end
if (nargout < 9)
    disp 'Note: skipping vector distribution. Communication statistics will display -1.';
    [I, s, p, q, r, c, rh, ch] = MatlabMondriaan(Pass, NumProcessors, Imbalance, -1, Permutation, Symmetry>0,SplitStrategy);
else
    if (nargout < 10)
        [I, s, p, q, r, c, rh, ch, B] = MatlabMondriaan(Pass, NumProcessors, Imbalance, -1, Permutation, Symmetry>0, SplitStrategy);
    else
        [I, s, p, q, r, c, rh, ch, B, u, v] = MatlabMondriaan(Pass, NumProcessors, Imbalance, -1, Permutation, Symmetry>0,SplitStrategy);
    end
    if (Symmetry == 2)
        B=A(p,q);
    end
end

end

