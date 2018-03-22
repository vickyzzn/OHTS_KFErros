
for b = 1:prog_pts
    
    for n = 1:N
        for v = 1:V
            for m = 1:3
               JK_error_OHTS_pre(n,v,b,m) = squeeze(JK_error_OHTS(n,v,OHTS_prog_loc{b}-1,m)); %Obtain MSE for each CV fold
                JK_error_OHTS_prog(n,v,b,m) = squeeze(JK_error_OHTS(n,v,ids(b),m));
            end            
        end
    end
    
end

