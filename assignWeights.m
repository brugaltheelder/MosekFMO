fid = fopen('singleMaskValue.bin','r');
singleMV=fread(fid,nVox,'int32');
fclose(fid);
% also load structSize


%%
over = zeros(1,29);
under= zeros(1,29);
threshHolder = zeros(1,29);
%%
structSize = zeros(1,29);
for s=1:29
    structSize(s) = length(find(log2(singleMV)==s-1));
    if structSize(s) == 0
        structSize(s) = 1
    end
    
end


%%
for i = 1:nVox
   aOver(i) = over(int32(log2(singleMV(i)))+1)/structSize(int32(log2(singleMV(i)))+1);
   aUnder(i) = under(int32(log2(singleMV(i)))+1)/structSize(int32(log2(singleMV(i)))+1);
   thresh(i) = threshHolder(int32(log2(singleMV(i)))+1);    
end

