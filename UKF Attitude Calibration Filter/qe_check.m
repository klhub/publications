qe_norm = zeros(1,tf/dts);

wait = waitbar(0,'Checking estimated quaternion ...');
for i=1:tf/dts-1
    qe_norm(:,i)=norm(qe(:,i));    
    waitbar(i/(tf/dts),wait)    
end
close(wait)
