function [start,fin] = check(data,offset)

flag = 0;
go = 0;
if(offset < 40)
    offset = offset + 40;
end
start = offset;
value = data(offset);
for i=offset:length(data)-40
    if(flag == 0 && go == 0)
        for j=-40:40
            if(data(i+j)>value+0.4 || data(i+j)<value-0.4 || data(i)<0.3)
                flag =1;
                break;
            else
                start = i-40;
                k = i+40;
                while(data(k)<value+0.8 && data(k)>value-0.8 && k<length(data))
                    k = k+1;
                end
                fin = k;
                go = 1;
                break;
            end
        end
    elseif(flag == 1 && go == 0)
        value = data(i);
        start = i;
        flag = 0;
    end
end