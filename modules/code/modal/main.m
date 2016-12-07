function main()
%% main
% 
% 
% 
% author: john devitis
% create date: 29-Nov-2016 09:48:39
	
	load grid1_nfa
    
    % extract vertical modal dof and natural frequencies
    V = squeeze(nfa.U(:,3,:));
    F = nfa.freq;
    
    % get FRF matrix of size [no x ni x ns]
    %  note - without specifiying input & output dof, default to all
    % HH - frf
    % AA - 
    [HH,AA,Wn,Qr,root,w] = nfa2frf(V,F); 
	
	% plot FRF measurement with contribution from each mode
    in = 6; out = 6;
    [hs,hh] = vibsFRF(AA,root,in,out,w);
    fh = vibsFRFplot(hs,hh,in,out,w);
    
    % plot equivalent measurement in time domain
    %  this is the impulse response function formed form the same residue
    %  matrix
    fs = 200; % sampling frequency [hz]
    l = 1;    % length of time vector (seconds)
    [hs,h] = vibsIRF(AA,root,in,out,fs,l);
    fh = vibsIRFplot(hs,h,in,out,fs,l);
    
    %% save FRF to file
    % reshape 3D array to 2D
    numcolumns = size(HH,1)*size(HH,2);
    FRF_R = reshape(permute(real(HH),[3 1 2]),[],numcolumns);
    FRF_I = reshape(permute(imag(HH),[3 1 2]),[],numcolumns);
    THead = 'Out%dIn%d%s';
    RI = ['R'; 'I'];
    n=1;
    for kk = 1:2
        for ii =  1:size(HH,2)
            for jj = 1:size(HH,1)
                labels{n} = sprintf(THead,[jj,ii,RI(kk)]);
                n=n+1;
            end
        end
    end
    
    % add ns values
    labels = [{'spectral'} labels];
    % write to table
    FRF_t = array2table(horzcat(w', FRF_R, FRF_I),'VariableNames',labels);
    % write table to text file
    writetable(FRF_t,'mod1_FRF.csv','Delimiter',',')
    
    dlmwrite('FRF_mod1',
end
