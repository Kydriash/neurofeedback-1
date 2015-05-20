function w = calc_csp_filter(eegdata1,eegdata2)
%%%%% eegdata1(n,t1), eegdata2(n,t2) ---> w'
%%%%% n - number of signals, t1,t2 - numbers of samples

sp_covariance1 = eegdata1*eegdata1'/(sum(diag(eegdata1*eegdata1')));
sp_covariance2 = eegdata2*eegdata2'/(sum(diag(eegdata2*eegdata2')));
composite_sp_covariance = sp_covariance1+sp_covariance2; %%matrix_of_eigenvectors*diag_matrix_of_eigenvalues*matrix_of_eigenvectors'
[eigenvectors, eigenvalues] = eig(composite_sp_covariance);
%sort values in desc order
[eigenvalues, indices] = sort(diag(eigenvalues),'descend');
eigenvectors = eigenvectors(:,indices);
whitening_transf_mx = sqrt(pinv(diag(eigenvalues)))*eigenvectors';
av_covariance_mx1 = whitening_transf_mx * sp_covariance1 * whitening_transf_mx';
av_covariance_mx2 = whitening_transf_mx * sp_covariance2 * whitening_transf_mx';

%%%check if sum of eigenvalues of av_covariance_mx1 and av_covariance_mx2
%%%== eye
[U1,Psi1] = eig(av_covariance_mx1); %#ok<ASGLU>
[U2,Psi2] = eig(av_covariance_mx2); %#ok<ASGLU>
disp(Psi1+Psi2)% == eye(length(Psi1)));
[B, D] = eig(av_covariance_mx1,av_covariance_mx2);

[D, indices] = sort(diag(D)); %#ok<ASGLU>
B = B(:,indices);
w = B'*whitening_transf_mx;

end