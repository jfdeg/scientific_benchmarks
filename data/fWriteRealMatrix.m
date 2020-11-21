function fWriteRealMatrix(sig_in,path,N)

        % real
        binFile = fopen([path num2str(N) '.bin'], 'w' );
        fwrite( binFile, reshape(real(sig_in), numel(real(sig_in)),1), 'double' );
        fclose( binFile );

end