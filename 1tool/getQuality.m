function Q = getQuality(path)
    jobj=jpeg_read(path);    
    q_table=jobj.quant_tables{1};
    
    % º∆À„Quality
    Q100 = q_table(end,6);
    if(Q100==1)
        Q = 100;
    elseif(Q100 < 100)
        Q = ( 2 - (Q100*0.01) ) * 50;
    else
        Q = (50*100)/Q100;
    end
end