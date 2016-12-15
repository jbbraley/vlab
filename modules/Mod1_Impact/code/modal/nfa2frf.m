function [HH,AA,Wn,Qr,root,w] = nfa2frf(V,F,in,out)
%% nfa2frf
% 
% 
% 
% author: john devitis
% create date: 29-Nov-2016 09:09:35

    %% log events
    logg = logger('NFA->FRF');
    logg.print('Setting things up...')
    
    %% Setup 
    if nargin < 3; in = 1:size(V,1);  end
    if nargin < 4; out = 1:size(V,1); end
    
    ni = length(in);    % number of inputs
    no = length(out);   % number of outputs
    ne = size(V,2);     % number of effective modes
	
    % radial natural frequency
    W = 2*pi*F; % [rad/sec]
    
    % form modal mass and stiffness
    %  assuming nfa produces true (mass scaled) eigenvectors. this produces
    %  a modal mass matrix of unity.
    Mr = eye(ne);                   % modal mass
    Kr = bsxfun(@times,eye(ne),W);  % modal stiffness
    
    % add proportional damping - 5% critical
    %  assuming damped natural freq roughly equal to undamped natural freq
    percd = .05;
    dampr = ones(length(F),1)*percd; % damping ratio [% critical per mode]
    dampf = -dampr.*W;               % damping factor [rad/sec]
    Wn = sqrt(W.^2 + dampf.^2);      % damped natural frequency [rad/sec]
    root = dampf + 1j.*Wn;           % form positive poles 
    
    % define modal scaling for unit mass
    %  notes: -true unity mass due to mass normalize eigenvector 
    %         -5 modes solved for
    Qr = 1./(2j.*diag(Mr).*Wn);
	
    % form radial frequency vector
    ns = 2^9;                   % number of spectral samples
    w = linspace(-200,200,ns);  % [rad/sec]
    
	%% Form FRF via residues
    logg.print('Forming FRF from Residues');
    
    AA = zeros(no,ni,ns); 
    HH = zeros(no,ni,ns);
    
    for ii = 1:ne   % loop modes    
        
        logg.task('Mode',ii,ne)
        
        % form [A] for mode ii -> [no x ni x ns]
        AA(:,:,ii) = Qr(ii) * V(:,ii) * V(:,ii)';     
        for jj = 1:no       % loop outputs
            for kk = 1:ni   % loop inputs            
                outid = out(jj); % output DOF index
                inid  = in(kk);  % input  DOF index            
                for ll = 1:ns % loop spectral lines         

                    % form [H] - add mode ii contribution -> [no x ni x ns]
                    %            complex conjugate
                    tt = AA(outid,inid,ii) ./ (1j*w(ll) - root(ii)) + ...
                              conj(AA(outid,inid,ii))./(1j*w(ll) - conj(root(ii)));

                    % add each mode for total response
                    HH(jj,kk,ll) = HH(jj,kk,ll) + tt;  
                end            
            end
        end   
        
        logg.done()
    end
    
    logg.finish()

end % /nfa2frf

