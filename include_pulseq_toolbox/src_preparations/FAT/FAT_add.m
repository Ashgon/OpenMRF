%% add fat suppression
if strcmp(FAT.mode, 'on')
    seq.addBlock(FAT.rf);
    seq.addBlock(FAT.gz_crush);
end