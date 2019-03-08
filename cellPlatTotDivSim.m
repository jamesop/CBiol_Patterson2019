function [cellTotPar,matCellGrowth]=cellPlatTotDivSim(growthRate,numCells,simTime,results,uDivLens,gDivSize)

threshTime=0;

guessMean=gDivSize;
matCellGrowth=normrnd(guessMean,abs(0.1*guessMean),[numCells,1]);


for i=1:simTime;
    i
    indMax=sum(any(~isnan(matCellGrowth),1));
    matCellGrowth=matCellGrowth(:,1:indMax);
    matGrow=~any(matCellGrowth==0,2);
    %find the last valid point of a cell
    currCellEnd=sum(~isnan(matCellGrowth),2);
    %check which of these last points are valid (is cell still growing)
    valEnds=currCellEnd(matGrow);
    
    %valid growing cell names;
    valGrowCells=find(matGrow);
    %find growth vals that are valid
    indGrowEnd=sub2ind(size(matCellGrowth),valGrowCells,valEnds);
    indNewGrowEnd=sub2ind([size(matCellGrowth,1),size(matCellGrowth,2)+1],valGrowCells,valEnds+1);
    
    %Current growing or dividing ends
    currSizes=matCellGrowth(indGrowEnd);
    
    
    %Generate division probabilities
    currpDivEv=hillFunc(results,currSizes);
    %Include if using a timer
    currpDivEv=currpDivEv.*(valEnds>threshTime).*(currSizes>min(uDivLens));
    %Generate growth amounts
    currGrow=currSizes+currSizes.*(normrnd(growthRate,2.2384e-04)); %0.006 works well
    
    %generate random numbers to choose which cells grow/don't
    matRands=rand(numel(currSizes),1);
    currDiv=matRands<currpDivEv;
    %note a division, by naming it 0
    currGrow(currDiv)=0;
    
    %make new horz cat to deal with current cells:
    allCurrCells=nan(size(matCellGrowth,1),1);
    allCurrCells(logical(matGrow))=currGrow;
    matCellGrowth=[matCellGrowth,nan(size(matCellGrowth,1),1)];
    matCellGrowth(indNewGrowEnd)=currGrow;
    
    %count number of division, mult by 2 to initiate new cells
    numDiv=sum(currDiv);
    
    matNewPar=[];
    if numDiv>0
        valPars=valGrowCells(currDiv);    
        matNewPar=nan(2*numDiv,2);
        %double up cells
        matNewPar(1:2:end,1)=valPars;
        matNewPar(2:2:end,1)=valPars;
        newChildNames=[size(matCellGrowth,1)+1:size(matCellGrowth,1)+numDiv*2];

        %new Children
        matNewPar(1:end,2)=newChildNames';

        %new cell holder
        newCurrCells=nan(numDiv*2,size(matCellGrowth,2));

        %calcBirthSize
        newBirthSizes=nan(numDiv*2,1);
        baseSize=(currSizes(currDiv)/2)*1.1;
        assym=normrnd(0.5,0.1);
        newBirthSizes(1:2:end)=2*baseSize*assym;
        newBirthSizes(2:2:end)=2*baseSize*(1-assym);
        
        newCurrCells(:,1)=newBirthSizes;
        matCellGrowth=[matCellGrowth;newCurrCells];
    end
    matTotPar=[valGrowCells(~currDiv),valGrowCells(~currDiv)];
    matTotPar=[matTotPar;matNewPar];
    
    cellTotPar{i}=matTotPar;
    
end


function outVals=hillFunc(inVals,x)

min=inVals(1);
max=inVals(2);
ec50=inVals(3);
hillc=inVals(4);

x1=x(:,1);

outVals=min+(max-min)./(1+(x1/ec50).^hillc);

