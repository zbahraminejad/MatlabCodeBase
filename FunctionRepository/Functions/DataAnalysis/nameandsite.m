function wellnameandsite = nameandsite(shot)
location = shot;
row = str2num(location(1));
if length(location)==6
    col = str2num([location(3) location(4)]);
else
    col = str2num(location(3));
end
site = str2num(location(end));




if row == 1
    rowletter = 'A';
elseif row == 2
    rowletter = 'B';
elseif row  == 3
    rowletter = 'C';
elseif row  == 4
    rowletter = 'D';
elseif row  == 5
    rowletter = 'E';
elseif row  == 6
    rowletter = 'F';
elseif row  == 7
    rowletter = 'G';
elseif row  == 8
    rowletter = 'H';
end

if col<10
    colnum = ['0' num2str(col)];
else
    colnum = num2str(col);
end 

wellnameandsite = [rowletter, colnum,'\site_',num2str(site),'\'];

end