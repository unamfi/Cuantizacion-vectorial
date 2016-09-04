%% K-MEDIAS

function centroides = kMedias(vectores,k)

    for i = 1:k
        z{i} = vectores{i};
        zAnt{i} = zeros(1,length(vectores{i}));
    end
    
    while notEqual(z,zAnt,k)
        
        zAnt = z;
        
        display('entro')
        
        for vector = 1:length(vectores)
            for cent = 1:length(z)
                d{vector}(cent) = sum((z{cent}-vectores{vector}).^2);
            end
        end
        
        vect=cell(1,k);
        for vector = 1:length(vectores)
            
            [val,ind]=min(d{vector});
            display(ind)
            temp=[vectores{vector};vect{ind}];
            
            vect{ind}=temp;
                
        end
        
        for n = 1:k
            z{n} = mean(vect{n},1);
        end
        
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