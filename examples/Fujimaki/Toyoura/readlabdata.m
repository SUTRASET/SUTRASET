    fn=fopen('LAB.DAT');
    sno=str2double(fgetl(fn));
    cno=str2double(fgetl(fn));
    slab=zeros(2,30);
    clab=zeros(2,30,sno);
% salinity

    for i=1:3
        for j=1:5
         ln(i,j)={fgetl(fn)};
        end
        temp=str2double(ln(i,4));
        slab(1:2,1:temp,i)=fscanf(fn, '%g %g', [2 temp]);
    end
% concentration
    for i=1:3
        for j=1:5
         ln2(i,j)={fgetl(fn)};
        end
        temp=str2double(ln2(i,4));
        clab(1:2,1:temp,i)=fscanf(fn, '%g %g', [2 temp]);
    end
% surface evaporation
    for i=1:5
    ln(i)={fgetl(fn)};
    end
    temp=str2double(ln(4));
    eslab=fscanf(fn, '%g %g', [2 temp]);
% ground evaporation
    for i=1:5
    ln(i)={fgetl(fn)};
    end
    temp=str2double(ln(4));
    eglab=fscanf(fn, '%g %g', [2 temp]);   
fprintf(1,'Lab.dat reading finished\n');
