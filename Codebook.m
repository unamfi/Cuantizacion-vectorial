%% CodeBook

function LPCCTotales = Codebook(s, figura)
    
    S1 = [];
    S2 = [];
    S3 = [];
    S4 = [];

    for i=1:length(s)
        display(s{i})
        LPCCres=LPCC(s{i});
        %Dividir pronunciación en 4 partes
        tamaniodesegmento = int8(length(LPCCres)/4);
        for j=1:tamaniodesegmento - 1
           S1 = [S1;LPCCres{j}];
        end
        for j=tamaniodesegmento:2 * tamaniodesegmento - 1 
           S2 = [S2;LPCCres{j}];
        end
        for j=2 * tamaniodesegmento:3 * tamaniodesegmento - 1 
           S3 = [S3;LPCCres{j}];
        end
        for j=3 * tamaniodesegmento:length(LPCCres) 
           S4 = [S4;LPCCres{j}];
        end
    end
    
    x = size(S1);
    
    for i=1:x(1)
        S1Cell{i} = S1(i,:);
    end
    
    x = size(S2);
    
    for i=1:x(1)
        S2Cell{i} = S2(i,:);
    end
    
    x = size(S3);
    
    for i=1:x(1)
        S3Cell{i} = S3(i,:);
    end
    
    x = size(S4)
    
    for i=1:x(1)
        S4Cell{i} = S4(i,:);
    end
    
%     for i=1:length(S1)
%         S1Cell{i} = S1(i,:);
%     end
%     
%     for i=1:length(S2)
%         S2Cell{i} = S2(i,:);
%     end
%     
%     for i=1:length(S3)
%         S3Cell{i} = S3(i,:);
%     end
%     
%     for i=1:length(S4)
%         S4Cell{i} = S4(i,:);
%     end

    centroides = kMedias(S1Cell, 16);
    save(strcat('centroides',figura,'1'),'centroides');
    centroides = kMedias(S2Cell, 16);
    save(strcat('centroides',figura,'2'),'centroides');
    centroides = kMedias(S3Cell, 16);
    save(strcat('centroides',figura,'3'),'centroides');
    centroides = kMedias(S4Cell, 16);
    save(strcat('centroides',figura,'4'),'centroides');
    

   
    %display(LPCCTotales)

end

%% LPCC

function LPCCres = LPCC(s)
    [R LPC p] = lpc(s);
    Q = 3*p/2;
    for i=1:1:length(LPC)
       
        C = zeros(1,Q+1);
        C(1) = log(sqrt(R{i}(1)));
        
        for m=1:Q
            if m>=1 && m<=p
                C(m+1)=LPC{i}(m);
                for k=1:m-1
                    C(m+1)=C(m+1)+((k/m)*C(k+1)*LPC{i}(m-k));
                end
            end
            if m>p
                C(m+1)=0;
                for k=m-p:m-1
                    C(m+1)=C(m+1)+((k/m)*C(k+1)*LPC{i}(m-k));
                end
                
            end
        end
        
        LPCCres{i} = C;
    end

end

function [R LPC p] = lpc(s)
% LPC Example of a local function.
    fid = fopen(s, 'r');
    senal0 = fread(fid, inf, 'int16');
    %senal0 = Recortar(senal0);%%%%%%%%%%
    b = [1 -0.95];
    a = [1];
    p = 8;
    senal0 = filter(b, a, senal0);
    
    lenX = length(senal0);
    numMuestrasVentana = 128;
    nV = lenX/numMuestrasVentana;
    
    %display(lenX)
    
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

end

%% K-MEDIAS

function centroides = kMedias(vectores,k)

    for i = 1:k
        z{i} = vectores{i};
        zAnt{i} = zeros(1,length(vectores{i}));
    end
    
    while notEqual(z,zAnt,k) 
        
        zAnt = z;
        
        %display('entro')
        
        for vector = 1:length(vectores)
            for cent = 1:length(z)
                %fprintf('z')
                %display(z{cent})
                %fprintf('vector')
                %display(vectores{vector})
                d{vector}(cent) = sum((z{cent}-vectores{vector}).^2);
            end
        end
        
        vect=cell(1,k);
        for vector = 1:length(vectores)
            
            [val,ind]=min(d{vector});
            %display(ind)
            temp=[vectores{vector};vect{ind}];
            
            vect{ind}=temp;
                
        end
        
        for n = 1:k
            
            display(vect{n})    
            z{n} = mean(vect{n},1);
            if isempty(z{n})
                z = zAnt;
                break
            end
        end
        
        display('caca')
        display(z)
        
        centroides = z;
        %display(d)
        
    end
    
    centroides = z;
    %display(z)
    %centroides=0;

end

function res = notEqual(z,zAnt,k)
    
    res = 1;

    for s = 1:k
        b{s} = isequal(z{s},zAnt{s});
        res = res & b{s};
    end
    
    if res == 0
        res = 1;
    else 
        res = 0;
    end        
end