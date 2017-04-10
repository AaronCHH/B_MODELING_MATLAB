trnData = [X Y];
numMFs = 2;
mfType = 'gbellmf';
in_fis = genfis1(trnData,numMFs,mfType);
out_fis = anfis(trnData,in_fis);
answer=evalfis(X,out_fis)
