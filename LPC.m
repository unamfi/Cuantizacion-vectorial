function [R LPC p] = LPC(s)

    fid = fopen(s, 'r');
    senal0 = fread(fid, inf, 'int16');
    %senal0 = Recortar(senal0);
    b = [1 -0.95];
    a = [1];
    p = 8;
    senal0 = filter(b, a, senal0);
    
    lenX = length(senal0);
    numMuestrasVentana = 128;
    nV = lenX/numMuestrasVentana;
    
    display(lenX)
    
    for trama=1:nV

        %senal = senal(1:128);
        senal = senal0((trama*128)-127:(trama*128));
        senal = senal.*hamming(128);
        r = zeros(p+1,1);
        for k=0:p
            for m=1:127-k+1
                r(k+1) = r(k+1) + senal(m)*senal(m+k);
            end
        end

        R{trama} = r;

        %display(r); %Hasta aqui esta bien.

        e = zeros(p+1,1);
        e(1,1) = r(1);
        k = zeros(p,1);
        alpha = zeros(p,p);
        for i=1:p
            suma = 0;
            for j=1:(i-1) 
                suma = suma + alpha(j,i-1)*r(abs(i-j+1));
            end

            k(i) = (r(i+1) - suma) / e(i,1);
            %fprintf('k(%i) = %i\n', i, k(i) )

            alpha(i,i) = k(i);
            %display(alpha(i,i))
            for j=1:i-1
                alpha(j,i)= alpha(j,i-1) - k(i)*alpha(i-j,i-1);
                %display(alpha)
            end
            e(i+1,1) = (1-k(i)*k(i))*e(i,1); %esta bien
            %fprintf('e(%i) = %i\n', i + 1, k(i) )
        end

        %display(alpha)

        am = zeros(p,1);
        for m=1:p
            am(m) = alpha(m,p);
        end

        LPC{trama} = am;

        %display(am)

    end





