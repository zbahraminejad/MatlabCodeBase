function testarea(bwimage)
    
    cellareas = regionprops(bwimage,'Area');
    cellareas = [cellareas.Area];
    hist(cellareas,50); 


end
