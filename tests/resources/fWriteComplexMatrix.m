function fWriteComplexMatrix(sig_in,path,N)

sig_in = reshape(sig_in, numel(sig_in),1);
adjacent(1:2:2*length(sig_in)) = real(sig_in); %
adjacent(2:2:2*length(sig_in)+1) = imag(sig_in);
binFile = fopen([path num2str(N) '.bin'], 'w' );
fwrite( binFile, adjacent, 'double' );
fclose( binFile );

end