function nuc_mask=excludewarped(nuc_mask)
nuc_label=bwlabel(nuc_mask);
nuc_solidity=cell2mat(struct2cell(regionprops(nuc_label,'Solidity')));
warpedobjects=find(nuc_solidity<0.9);
nuc_label(ismember(nuc_label,warpedobjects))=0;
nuc_mask=nuc_label>0;