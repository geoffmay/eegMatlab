grow = grower;
a = [];

testSize = 1000000;

tic;
counter = 1;
for i = 1:testSize;
    aTime(counter) = 1;
    counter = counter + 1;
end
aTime = toc;

tic;
for i = 1:testSize;
grow.add(1);
end
gTime = toc;