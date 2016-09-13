function mondriaanplot(A, NumProcessors, Imbalance, Permutation, Symmetry, Iterations, SplitStrategy)
% mondriaanplot   Creates a coloured plot of the Mondriaan partitioning of
%                 a given matrix A over NumProcessors processors with a
%         maximum imbalance specified by Imbalance.
%         If desired, a permutation may also be specified in the
%         variable Permutation: 0 = no permutation,
%                               1 = move cut rows/columns to the front
%                               2 = move cut rows/columns to the middle
%                               3 = move cut rows/columns to the back.
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
%         It is also possible to set the desired maximum number of
%         Iterations to stop Mondriaan while still in recursing to make
%         the initial NumProcessors parts; this can be instructive as it
%         allows peeking in the intermediate steps Mondriaan takes.
%
%         This can also be used to make animations, but requires several
%         Mondriaan runs. Use of the MondriaanPlot command-line utility
%         should be prefered for that end.
%         Setting this value <= 0 will lead to a full splitting.
%         
%         Return value:
%             Nothing is returned, but a picture of the partitioning process is generated.
%
%         Usage:
%             mondriaanplot(A, NumProcessors, Imbalance, Permutation, Symmetry, Iterations, SplitStrategy);
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
    Symmetry = 0;
    Iterations = -1;
end
if (nargin < 5)
    Symmetry = 0;
    Iterations = -1;
end
if (nargin < 6)
    Iterations = -1;
end
if (nargin<7)
    SplitStrategy = -1;
elseif (SplitStrategy <-1 || SplitStrategy > 8)
    disp('Invalid value for SplitStrategy, falling back to the default one.')
    SplitStrategy = -1;
end

if( Symmetry < 0 || Symmetry > 2)
    disp('Invalid value for Symmetry');
    return;
elseif( Symmetry > 0 )
    disp('------------------------------------------------------------------------------------------------');
    disp('|Warning: the input matrix is assumed symmetric, this call will not fail if it actually is not!|');
    disp('------------------------------------------------------------------------------------------------');
    disp('Info: recommend symmetric split strategy. Other strategies are valid,');
    disp('      but may not yield what you may expect.');
    disp('Info: forcing symmetry overrides "UseSingleEntry"  option to "yes"');
    disp('Info: forcing symmetry overrides "SingleEntryType" option to "lower"');
end
if (Symmetry > 0)
    Pass = tril(A);
else
    Pass = A;
end

Colours = colormap(hsv(NumProcessors + 1));

[I, s, p, q, r, c, rh, ch, B, u, v] = MatlabMondriaan(Pass, NumProcessors, Imbalance, Iterations, Permutation, Symmetry>0, SplitStrategy);

% Create figure.
clf
hold off

for i=1:NumProcessors
    C = (I(p,q) == i);
    
    % The following code is based on spy.m, included with Matlab.
    %spy(C, Colours(i,:));
    [k, l] = find(C);
    plot(l, k, 'marker', '.', 'markersize', 2, 'linestyle', 'none', 'color', Colours(i,:));
    hold on
end

%Print some basic statistics.
[m, n] = size(A);
xlabel(['m = ', int2str(m), ', n = ', int2str(n), ', nz = ', int2str(nnz(A)), ', p = ', int2str(NumProcessors), ', eps = ', num2str(s(2)), ', max. com. = ', num2str(s(3)), ', com. vol. = ', num2str(s(4))]);

set(gca, 'xlim', [0 n+1], 'ylim', [0 m+1], 'grid', 'none', 'plotboxaspectratio', [n+1 m+1 1], 'ydir', 'reverse');

end
