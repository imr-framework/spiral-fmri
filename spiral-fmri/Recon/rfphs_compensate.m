 
function [corrected_kspace] = rfphs_compensate (k_space_dat)
%% Calculate the rf phase vector %%

rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 150;

for i=1:54*12
    rfphs_v(i) = rfphs; 
    rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
    rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
end

% move the first value to the last element of the vector
rfphs_v = circshift(rfphs_v, 1);
% receiver_phase = repmat([-pi pi], 1, 54*12/2);
%rfphs_v = rfphs_v+ receiver_phase;

%% Multiply the k-space data by the phase %%
count = 1;
for i=1:size(k_space_dat,3)
    for j=1:size(k_space_dat,2)
        % positive phase
        corrected_kspace(:,j,i,:)=k_space_dat(:,j,i,:)*exp(-1j*rfphs_v(count));
        count=count+1;
    end 
end

%% Align the spiral arms from different frames

rot = [1 2 3];

for fs_ksp = 1:12-2
    for sl = 1:54
        new_ksp(:,1,sl,fs_ksp,:) = corrected_kspace(:,sl,(fs_ksp-1)+find(rot==1),:);
        new_ksp(:,2,sl,fs_ksp,:) = corrected_kspace(:,sl,(fs_ksp-1)+find(rot==2),:);
        new_ksp(:,3,sl,fs_ksp,:) = corrected_kspace(:,sl,(fs_ksp-1)+find(rot==3),:);
        rot = circshift(rot,-1);
    end
    
end

corrected_kspace = new_ksp;



end
