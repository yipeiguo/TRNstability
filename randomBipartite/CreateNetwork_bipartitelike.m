% This is a function to create a bipartite-like regulatory network by
% specfying the number of TFs, the total number of regulatory interactions,
% and the fraction of which that are TF-nonTF.

function IntTypeVec = CreateNetwork_bipartitelike(NetworkProp)

    % Extract network parameters
    N = NetworkProp.N;
    numTFs = NetworkProp.numTFs; % numnonTFs = N-numTFs;
    numInts = NetworkProp.numInts;
    numSelfInts = NetworkProp.numSelfInts; numNSints = numInts - numSelfInts;
    numInts_TFTF = NetworkProp.numInts_TFTF;
    numInts_TFnonTF = numNSints - numInts_TFTF;
%     numInts_TFnonTF = NetworkProp.numInts_TFnonTF;
%     numInts_TFTF = numInts - numInts_TFnonTF;
    if isfield(NetworkProp,'Pneg_selfint')
        Pneg_selfint = NetworkProp.Pneg_selfint;
    else
        Pneg_selfint = 0.5;
    end
    if isfield(NetworkProp,'Pneg_nonselfint')
        Pneg_nonselfint = NetworkProp.Pneg_nonselfint;
    else
        Pneg_nonselfint = 0.5;
    end
    
    
    TFinds = (1:numTFs)';
    nonTFinds = (numTFs+1:N)';
    
    % Draw self interactions:
    selfregInds = randperm(numTFs,numSelfInts)';        
    
    % Draw non-self interactions:
    IntPair_TFinds = randsample(TFinds,numNSints,true);
    % Draw TF-TF interactions:
    IntPair_regulatedinds_TF = randsample(TFinds,numInts_TFTF,true);
    IntPairs_TFTF = [IntPair_TFinds(1:numInts_TFTF),IntPair_regulatedinds_TF];
    IntPair_nonself = IntPairs_TFTF(IntPairs_TFTF(:,1)~=IntPairs_TFTF(:,2),:);
    numnonselfdrawn = size(IntPair_nonself,1);
    [~,ia] = unique(IntPair_nonself,'rows');
    while length(ia) < numInts_TFTF
        if numnonselfdrawn < numInts_TFTF
            numtoadd = numInts_TFTF-numnonselfdrawn;
            IntPair_toadd = [randsample(TFinds,numtoadd,true),...
                randsample(TFinds,numtoadd,true)];
        else
            IntPair_toadd = [];
        end
        repeatInd = setdiff(1:numnonselfdrawn,ia);
        IntPair_nonself(repeatInd,:) = [randsample(TFinds,length(repeatInd),true),...
            randsample(TFinds,length(repeatInd),true)];
        IntPairs_TFTF = [IntPair_nonself;IntPair_toadd];
            
        IntPair_nonself = IntPairs_TFTF(IntPairs_TFTF(:,1)~=IntPairs_TFTF(:,2),:);
        numnonselfdrawn = size(IntPair_nonself,1);
        [~,ia] = unique(IntPair_nonself,'rows');
    end           
    
    % Draw TF-nonTF interactions:
    IntPair_regulatedinds_nonTF = randsample(nonTFinds,numInts_TFnonTF,true);
    IntPairs_TFnonTF = [IntPair_TFinds(1+numInts_TFTF:end),IntPair_regulatedinds_nonTF];
    [~,ia] = unique(IntPairs_TFnonTF,'rows');
    while length(ia) < numInts_TFnonTF
        repeatInd = setdiff(1:numInts_TFnonTF,ia);
        IntPairs_TFnonTF(repeatInd,:) = [randsample(TFinds,length(repeatInd),true),...
            randsample(nonTFinds,length(repeatInd),true)];
        [~,ia] = unique(IntPairs_TFnonTF,'rows');
    end           
    % combine all interaction pairs:
    IntPairs_self = [selfregInds,selfregInds];
    IntPairs_all = [IntPairs_self;IntPairs_TFTF;IntPairs_TFnonTF];
    if length(unique(IntPairs_all(:,1))) < numTFs
        warning('numTFs less than specified!');
    end
    
    % specify interaction type:
    IntInds = sub2ind([N,N],IntPairs_all(:,1),IntPairs_all(:,2));
    IntTypeVec = zeros(N*N,1);
    IntTypeVec(IntInds(1:numSelfInts)) = sign(rand(numSelfInts,1)-Pneg_selfint);
    IntTypeVec(IntInds(1+numSelfInts:end)) = sign(rand(numInts-numSelfInts,1)-Pneg_nonselfint);

%     % Draw TF-TF interactions;
%     Tmat = sign(sprandn(numTFs,numTFs,numInts_TFTF/(numTFs^2)));
%     
%     % Draw TF-nonTF interactions:
%     Qmat = sign(sprandn(numTFs,numnonTFs,numInts_TFnonTF/(numTFs*numnonTFs)));
%     
%     IntTypeMat = [[Tmat,Qmat];zeros(numnonTFs,N)];
%     IntTypeVec = IntTypeMat(:);
        
end