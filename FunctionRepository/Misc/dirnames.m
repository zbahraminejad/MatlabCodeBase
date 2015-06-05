function filenames = dirnames(directory)
names = dir(directory);
names = {names.name};
names(1:2) = [];
filenames = names;
end