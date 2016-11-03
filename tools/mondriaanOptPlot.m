function mondriaanOptPlot(A, Imbalance, Volume)
% mondriaanplot   Creates a coloured plot of the MondriaanOpt partitioning of
%                 a given matrix A over 2 processors. The computed
%                 partitioning has the lowest communication volume among
%                 all partitionings that have imbalance at most epsilon.
%                 
%         Required arguments:
%             A         = Sparse matrix
%             Imbalance = Maximum allowed imbalance
%             Volume    = The initial upper bound for the communication volume
%
%         Return value:
%             Nothing is returned, but a picture of the partitioning process is generated.
%
%         Usage:
%             mondriaanOptPlot(A, Imbalance, Volume);
%
    if (nargin < 3)
        disp 'Error: the three arguments to mondriaanOptPlot function are mandatory!';
        return;
    end
    if (Imbalance < 0)
        disp('Invalid value for Imbalance');
        return;
    end
    if (issparse(A) ~= 1)
        disp('Input matrix A is not sparse');
        return;
    end

    NumProcessors = 2;

    Colours = [[1,0,0];[0,0,1];[1,1,0]];

    [I, s] = mondriaanOpt(A, Imbalance, Volume);

    % Create figure.
    clf
    hold off

    for i=1:NumProcessors+1
        C = (I == i);

        [k, l] = find(C);
        for t=1:length(k)
            rectangle('Position',[l(t)-0.5 k(t)-0.5 1 1], 'FaceColor', Colours(i,:), 'EdgeColor', 'none');
        end
        hold on
    end
    
    % Draw vectors
    [m,n] = size(A);
    for i=1:m
        c = (sum(I(i,:)==1)>0)*1 + (sum(I(i,:)==2)>0)*2;
        if(sum(I(i,:)==3)>0)
            c = 3; % Force c to be three.
            % It can happen that I(i,.)==3 while there is only one
            % processor represented in the row, and all other nonzeros are
            % free. This happens whenever MondriaanOpt cannot distribute
            % these free nonzeros to a single processor while also staying
            % within load imbalance limits.
        end
        if c > 0
            rectangle('Position',[0 i-0.5 0.5 1], 'FaceColor', Colours(c,:));
        end
    end
    for j=1:n
        c = (sum(I(:,j)==1)>0)*1 + (sum(I(:,j)==2)>0)*2;
        if(sum(I(:,j)==3)>0)
            c = 3; % Force c to be three. (See above)
        end
        if c > 0
            rectangle('Position',[j-0.5 0 1 0.5], 'FaceColor', Colours(c,:));
        end
    end

    %Print some basic statistics.
    [m, n] = size(A);
    xlabel(['m = ', int2str(m), ', n = ', int2str(n), ', nz = ', int2str(nnz(A)), ', p = ', int2str(NumProcessors), ', eps = ', num2str(s(2)), ', com. vol. = ', num2str(s(3))]);

    set(gca, 'xlim', [0 n+1], 'ylim', [0 m+1], 'plotboxaspectratio', [n+1 m+1 1], 'ydir', 'reverse');

end
