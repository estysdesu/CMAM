function [ octreeList] = uniformOctreeDepth(octreeCells, octreeLength)
    % Description:
    %   uniformOctreeDepth creates a uniform list of octree sequences that have the same depth. Full nodes that terminate higher in the tree than the lowest level, are expanded into the lowest branch of the tree with all their octree node components full as well.
    % Inputs:
    %   octreeList: List of octree sequences having the format r<0-9>*
    % Outputs:
    %   newOctreeLIst: List of octree sequences having the same format as the input, but ensured that all the octree sequences will be of equal length

    % calculate max length
    if nargin < 2
        MAX_LENGTH = 0;
        for n = 1:length(octreeCells)
            if length(octreeCells{n}) > MAX_LENGTH
            MAX_LENGTH = length(octreeCells{n});
            end
        end
    else
        MAX_LENGTH = octreeLength;
    end
    
    % convert cells to char arrays with nonuniform rows containing empty buffer chars (' ')
    octreeList = char(ones(length(octreeCells), MAX_LENGTH) * ' ');
    for n=1:length(octreeCells)
        m = length(octreeCells{n});
        octreeList(n, 1:m) = octreeCells{n};
    end
    
    % expand shorter leafs to full length of longest
    OCTREE_ELEMS = ['0'; '1'; '2'; '3'; '4'; '5'; '6'; '7'];
    while ~all(octreeList(:) ~= ' ')
        [r, c] = find(octreeList == ' ', 1);
        tmpOctreeList = char(ones(size(OCTREE_ELEMS)) * octreeList(r, 1:c-1)); % generate list of 8 identical r sequences
        tmpOctreeList(:, end+1) = OCTREE_ELEMS; % append new octree level
        
        padLen = MAX_LENGTH - length(tmpOctreeList(1, :));
        if padLen > 0
            tmpOctreeList = [tmpOctreeList, char(ones(size(tmpOctreeList, 1), padLen) * ' ')]; % append blanks
        end
            
        octreeList = [octreeList(1:r-1, :); tmpOctreeList; octreeList(r+1:end, :)];
    end

