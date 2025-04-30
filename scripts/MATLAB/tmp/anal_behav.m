subjlist = [1:3 5:13];
glm = 2;
for s = 1:length(subjlist)
    R = construct_dsgmat(sprintf('S%02d',subjlist(s)),glm);
    nRun = 8;
    if subjlist(s)~=7
        R2 = construct_dsgmat(sprintf('R%02d',subjlist(s)),glm);
        R = combine_behavdata(R,R2);
        nRun = 16;
    end
    for r=1:nRun
        temp = R.MT(r,R.isRepBoth(r,:)==1);
        MTr(s,r) = mean(temp(temp~=0));
        temp = R.MT(r,R.isRepBoth(r,:)==0);
        MTnb(s,r) = mean(temp(temp~=0));        
        temp = R.RT(r,R.isRepBoth(r,:)==1);
        RTr(s,r) = mean(temp(temp~=0));
        temp = R.RT(r,R.isRepBoth(r,:)==0);
        RTnb(s,r) = mean(temp(temp~=0));
        
        temp = R.MT(r,R.isRepBoth(r,:)==1 & R.isLetter(r,:));
        MTrl(s,r) = mean(temp(temp~=0));
        temp = R.MT(r,R.isRepBoth(r,:)==0 & R.isLetter(r,:));
        MTnbl(s,r) = mean(temp(temp~=0));        
        temp = R.RT(r,R.isRepBoth(r,:)==1 & R.isLetter(r,:));
        RTrl(s,r) = mean(temp(temp~=0));
        temp = R.RT(r,R.isRepBoth(r,:)==0 & R.isLetter(r,:));
        RTnbl(s,r) = mean(temp(temp~=0));        

        temp = R.MT(r,R.isRepBoth(r,:)==1 & R.isSpatial(r,:));
        MTrs(s,r) = mean(temp(temp~=0));
        temp = R.MT(r,R.isRepBoth(r,:)==0 & R.isSpatial(r,:));
        MTnbs(s,r) = mean(temp(temp~=0));        
        temp = R.RT(r,R.isRepBoth(r,:)==1 & R.isSpatial(r,:));
        RTrs(s,r) = mean(temp(temp~=0));
        temp = R.RT(r,R.isRepBoth(r,:)==0 & R.isSpatial(r,:));
        RTnbs(s,r) = mean(temp(temp~=0));    
        
        temp = R.MT(r,R.isLetter(r,:)==1);
        MTl(s,r) = mean(temp(temp~=0));
        temp = R.MT(r,R.isLetter(r,:)==0);
        MTv(s,r) = mean(temp(temp~=0));        
        temp = R.RT(r,R.isLetter(r,:)==1);
        RTl(s,r) = mean(temp(temp~=0));
        temp = R.RT(r,R.isLetter(r,:)==0);
        RTv(s,r) = mean(temp(temp~=0));
        

        temp = R.MT(r,R.isNRep(r,:)==1);
        MTn(s,r) = mean(temp(temp~=0));        
        temp = R.RT(r,R.isNRep(r,:)==1);
        RTn(s,r) = mean(temp(temp~=0));
        
    end 
    Prn_MT(s) = signrank(MTr(s,1:nRun)-MTn(s,1:nRun));
    Plv_MT(s) = signrank(MTl(s,1:nRun)-MTv(s,1:nRun));
    Prn_RT(s) = signrank(RTr(s,1:nRun)-RTn(s,1:nRun));
    Plv_RT(s) = signrank(RTl(s,1:nRun)-RTv(s,1:nRun));   
    Prn_MT2(s) = signrank(MTr(s,1:nRun)-MTnb(s,1:nRun));
    Prn_RT2(s) = signrank(RTr(s,1:nRun)-RTnb(s,1:nRun));
    
    Diffrn_MT(s) = mean(MTr(s,1:nRun)-MTn(s,1:nRun));
    Difflv_MT(s) = mean(MTl(s,1:nRun)-MTv(s,1:nRun));
    Diffrn_RT(s) = mean(RTr(s,1:nRun)-RTn(s,1:nRun));
    Difflv_RT(s) = mean(RTl(s,1:nRun)-RTv(s,1:nRun));
    Diffrn_MT2(s) = mean(MTr(s,1:nRun)-MTnb(s,1:nRun));
    Diffrn_RT2(s) = mean(RTr(s,1:nRun)-RTnb(s,1:nRun));

    mMTv(s) = mean(MTv(s,1:nRun));
    mMTl(s) = mean(MTl(s,1:nRun));
    mRTv(s) = mean(RTv(s,1:nRun));
    mRTl(s) = mean(RTl(s,1:nRun));
    
    mMTr(s) = mean(MTr(s,1:nRun));
    mMTn(s) = mean(MTn(s,1:nRun));
    mRTr(s) = mean(RTr(s,1:nRun));
    mRTn(s) = mean(RTn(s,1:nRun));
end


data = [mRTl mRTv mMTl mMTv]; 
figure;
subplot(2,2,1);
bar([mean(mMTv) mean(mMTl)]);
hold on;
errorbar([mean(mMTv) mean(mMTl)], [std(mMTv)/sqrt(12) std(mMTl)/sqrt(12)], 'k', 'linestyle', 'none');

