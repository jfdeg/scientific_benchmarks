
write = 0;

%% Covariance Matrix Estimation benchmark %%
disp('Covariance Matrix Estimation generation');

nb_run = 5;

for N=100:100:3000
    
    sig_in = complex(randn(N,N,nb_run),randn(N,N,nb_run));
    
    if write == 1
        % real
        fWriteRealMatrix(sig_in,'./CovMatEstimation/CME_Real_in_',N)
        % complex
        fWriteComplexMatrix(sig_in,'./CovMatEstimation/CME_Complex_in_',N)
    end
    
    covar_mat_r=zeros(N,N,nb_run); % real
    covar_mat_c=zeros(N,N,nb_run); % complex
    
    for i=1:nb_run
        M = squeeze(sig_in(:,:,i));
        covar_mat_r(:,:,i) = (1/N)*real(M)*real(M)';
        covar_mat_c(:,:,i) = (1/N)*M*M';
        clear M;
    end
    clear sig_in;
    
    % for writing
    if write == 1
        fWriteRealMatrix(covar_mat_r,'./CovMatEstimation/CME_Real_out_',N);
        fWriteComplexMatrix(covar_mat_c,'./CovMatEstimation/CME_Complex_out_',N);
    end
     
    clear covar_mat_r;
    clear covar_mat_c;
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

