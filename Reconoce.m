function clase=Reconoce(s)
    LPCCres=LPCC(s);
    S1 = [];
    S2 = [];
    S3 = [];
    S4 = [];
    
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
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%
    load('centroidesCirculo1')
    for i=1:length(S1Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S1Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS1 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%
    load('centroidesCirculo2')
    for i=1:length(S2Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S2Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS2 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%
    load('centroidesCirculo3')
    for i=1:length(S3Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S3Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS3 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     %%%%%%%%
    load('centroidesCirculo4')
    for i=1:length(S4Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S4Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS4 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%
    distanciaaCirculo = distanciaTotalS1 + distanciaTotalS2 + distanciaTotalS3 + distanciaTotalS4;
    
    %%%%%%
    %%%%%%
    
    %%%%%%%%
    load('centroidesTriangulo1')
    for i=1:length(S1Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S1Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS1 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%
    load('centroidesTriangulo2')
    for i=1:length(S2Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S2Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS2 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%
    load('centroidesTriangulo3')
    for i=1:length(S3Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S3Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS3 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     %%%%%%%%
    load('centroidesTriangulo4')
    for i=1:length(S4Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S4Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS4 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%
    distanciaaTriangulo = distanciaTotalS1 + distanciaTotalS2 + distanciaTotalS3 + distanciaTotalS4;
    
    %%%%%%
    %%%%%%

     %%%%%%%%
    load('centroidesCuadrado1')
    for i=1:length(S1Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S1Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS1 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%
    load('centroidesCuadrado2')
    for i=1:length(S2Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S2Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS2 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%
    load('centroidesCuadrado3')
    for i=1:length(S3Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S3Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS3 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     %%%%%%%%
    load('centroidesCuadrado4')
    for i=1:length(S4Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S4Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS4 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%
    distanciaaCuadrado = distanciaTotalS1 + distanciaTotalS2 + distanciaTotalS3 + distanciaTotalS4;
    
    %%%%%%
    %%%%%%

     %%%%%%%%
    load('centroidesEstrella1')
    for i=1:length(S1Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S1Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS1 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %%%%%%%%
    load('centroidesEstrella2')
    for i=1:length(S2Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S2Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS2 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%
    load('centroidesEstrella3')
    for i=1:length(S3Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S3Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS3 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     %%%%%%%%
    load('centroidesEstrella4')
    for i=1:length(S4Cell)
        for j=1:length(centroides)
            
            distanciatemporal = [distanciatemporal sum((S4Cell{i}-centroides{j}).^2)];
        end
        val = min(distanciatemporal);
        distancias = [distancias val];
        
    end
    distanciaTotalS4 = sum(distancias);
    clear distancias
    clear distanciatemporal
    
    %%%%%%%%%%%% Definicion distancias %%%%%%%%%%%%%
    
    distanciatemporal = [];
    distancias = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%
    distanciaaEstrella = distanciaTotalS1 + distanciaTotalS2 + distanciaTotalS3 + distanciaTotalS4;
    
    %%%%%%
    %%%%%%

    resultados = [distanciaaCirculo distanciaaEstrella distanciaaCuadrado distanciaaTriangulo];
    
    [val, ind] = min(resultados)
    
    if ind == 1 
        fprintf('Es un circulo')
        imshow('circulo.png')
    end
    if ind == 2 
        fprintf('Es una estrella')
        imshow('estrella.jpg')
    end
    if ind == 3 
        fprintf('Es un cuadrado')
        imshow('cuadrado.jpg')
    end
    if ind == 4 
        fprintf('Es un triangulo')
        imshow('triangulo.png')
    end
   
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

