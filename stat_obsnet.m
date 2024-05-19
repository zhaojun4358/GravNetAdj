function stat_obsnet
unObsJD(1,1) = jdObsN(1,1);
st3 = 0;
for i = 1:numJ
    if strcmp(unObsJD,jdObsN(i,1)) == zeros(length(unObsJD),1)
        st3 = st3 + 1;
        unObsJD(st3,1) = jdObsN (i,1);
    end
end
                   fprintf(fid,'number of instrument£º%d\n',length(yiObs));
                   for i = 1:length(yiObs)
                       fprintf(fid,'%s\n',unObsN{i,1});
                   end
                   fprintf(fid,'\n number of obsolute gravity point :%d\n',st3);
                   for i = 1:st3
                       fprintf(fid,'%s\n',unObsJD{i,1});
                    end
                   fprintf(fid,'\n total point:%d\n',length(unObsN));
                   for i = 1:length(unObsN)
                       fprintf(fid,'%s\n',unObsN{i});
                   end