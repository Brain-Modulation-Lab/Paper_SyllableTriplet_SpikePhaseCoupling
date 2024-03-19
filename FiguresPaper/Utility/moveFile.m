function moveFile(From, To)

try
    system(sprintf('mv %s %s', From,To));
catch
    system(sprintf('move %s %s', From,To));
end
end