function GrabarVoz()
    display('Comenzando a grabar voz.')
    display('Presiona cualquier tecla para comenzar a grabar la palabra "Círculo".')
    pause
    display('Di "Círculo" claramente.')
    for i=1:10
        Fs = 8000;
        y = wavrecord(5*Fs,Fs,'int16');
        x = double(y);
        %plot(x);
        %display(x);
        fidw = fopen(strcat('Circulo',int2str(i),'.raw'),'w');
        for n = i:length(x)
            fwrite(fidw,y(n),'int16',0,'l');
        end
        fclose(fidw);
        soundsc(x,11025)
        pause(5)
        display(strcat('Repetición',' ',int2str(i),' terminada.'))
        if i < 10
            display(strcat('Presiona cualquier tecla para la repetición',int2str(i+1)))
            pause
        else
            display('Fin de las repeticiones para la palabra "Círculo"')
        end        
    end
    
    display('Presiona cualquier tecla para comenzar a grabar la palabra "Triángulo".')
    pause
    display('Di "Triángulo" claramente.')
    for i=1:10
        Fs = 8000;
        y = wavrecord(5*Fs,Fs,'int16');
        x = double(y);
        fidw = fopen(strcat('Triangulo',int2str(i),'.raw'),'w');
        for n = i:length(x)
            fwrite(fidw,x(n),'int16',0,'l');
        end
        fclose(fidw);
        soundsc(x,11025)
        pause(5)
        display(strcat('Repetición',' ',int2str(i),' terminada.'))
        if i < 10
            display(strcat('Presiona cualquier tecla para la repetición',int2str(i+1)))
            pause
        else
            display('Fin de las repeticiones para la palabra "Triángulo"')
        end        
    end
    
    display('Presiona cualquier tecla para comenzar a grabar la palabra "Cuadrado".')
    pause
    display('Di "Cuadrado" claramente.')
    for i=1:10
        Fs = 8000;
        y = wavrecord(5*Fs,Fs,'int16');
        x = double(y);
        fidw = fopen(strcat('Cuadrado',int2str(i),'.raw'),'w');
        for n = i:length(x)
            fwrite(fidw,x(n),'int16',0,'l');
        end
        fclose(fidw);
        soundsc(x,11025)
        pause(5)
        display(strcat('Repetición',' ', int2str(i),' terminada.'))
        if i < 10
            display(strcat('Presiona cualquier tecla para la repetición',int2str(i+1)))
            pause
        else
            display('Fin de las repeticiones para la palabra "Cuadrado"')
        end        
    end
    
    display('Presiona cualquier tecla para comenzar a grabar la palabra "Estrella".')
    pause
    display('Di "Estrella" claramente.')
    for i=1:10
        Fs = 8000;
        y = wavrecord(5*Fs,Fs,'int16');
        x = double(y);
        fidw = fopen(strcat('Estrella',int2str(i),'.raw'),'w');
        for n = i:length(x)
            fwrite(fidw,x(n),'int16',0,'l');
        end
        fclose(fidw);
        soundsc(x,11025)
        pause(5)
        display(strcat('Repetición',' ', int2str(i),' terminada.'))
        if i < 10
            display(strcat('Presiona cualquier tecla para la repetición',int2str(i+1)))
            pause
        else
            display('Fin de las repeticiones para la palabra "Estrella"')
        end        
    end
    
    display('=======================================')
    display('         Fin de la grabación           ')
    display('=======================================')