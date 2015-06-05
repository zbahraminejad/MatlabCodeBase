function selectmodedisplay(edgemask)
global channel
dcm_obj=datacursormode(gcf);
datacursormode on;
response='on';
while response                                          %get next graph: 'enter'
    response=input('','s');
    switch response                             %view angie sensor: 'a','enter'
        case '1'
            channel=1;
        case '2'
            channel=2;
        case '3'
            channel=3;
        otherwise
            response=false;
            continue;
    end
    figure(gcf);
    delete(findall(gcf,'Type','hggroup'));          %must remove datapoints before changing program
    if edgemask
        set(dcm_obj,'Updatefcn',@datatip_image_channelx_2);
    else
        set(dcm_obj,'Updatefcn',@datatip_image_nomask);
    end
    %set(dcm_obj,'Updatefcn',@datatip_image_channelx_nude);
    while ~strcmp(response,'d')                     %exit image mode:'d','enter'
        response = input('','s');
    end
    currentfigurenum = gcf;
    currentfigurenum = currentfigurenum.Number;
    delete(findall(gcf,'Type','hggroup'));
    close(currentfigurenum+1);                                   %close the image window
    figure(gcf);
end