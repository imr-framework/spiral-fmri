rfphs = 0;  % rad
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 150;

for i=1:648
    rfphs_v(i) = rfphs; 
    rphsLast = rfphs;
    rfphs = rfphs + (rf_spoil_seed / 180 * pi) * rf_spoil_seed_cnt; % rad
    rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
end