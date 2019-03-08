function [cellTotPar,matCellGrowth]=cellSizerTotDivSim(growthRate,numCells,simTime,numNoise,meanDivLen)


threshTime=60;

guessMean=meanDivLen/2;
matCellGrowth=normrnd(guessMean,abs(0.1*guessMean),[numCells,1]);

matDivSize=zeros(numCells,1);
 
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
    %get dimensions right
    addZero=zeros(size(matCellGrowth,1)-size(matDivSize,1),1);
    matDivSize=[matDivSize;addZero];
    matDivSize(currCellEnd==1)=normrnd(meanDivLen,numNoise,sum(currCellEnd==1),1);
    cDivSize=matDivSize(matGrow);
    %Include if using a timer
    currpDivEv=cDivSize.*(valEnds>threshTime);
    %Generate growth amounts
    currGrow=currSizes+currSizes.*growthRate; %0.006 works well
    
    %generate random numbers to choose which cells grow/don't
    
    currDiv=currSizes>currpDivEv & currpDivEv~=0;
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
        newBirthSizes(1:2:end)=baseSize-normrnd(1,0.1);
        newBirthSizes(2:2:end)=baseSize+normrnd(1,0.1);
        
        newCurrCells(:,1)=newBirthSizes;
        matCellGrowth=[matCellGrowth;newCurrCells];
    end
    matTotPar=[valGrowCells(~currDiv),valGrowCells(~currDiv)];
    matTotPar=[matTotPar;matNewPar];
    
    cellTotPar{i}=matTotPar;
    
end


