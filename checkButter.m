%frequencies are defined by 

sampleRate = 128;
lowPass = 1;
highPass = 3;
order = 6;

lo = lowPass / (sampleRate/2);
hi = highPass / (sampleRate/2);
[a1 b1] = butter(order, [lo, hi]);
[a b c d] = butter(order, [lo, hi]);
[zero, pole, gain] = butter(order, [lo,hi]);

Hd = dfilt.df2t(b1,a1);

x = (1:200) .* (pi / 100);
slow = sin(x);
fast = sin(x.*10);
fastFilt = filter(Hd, fast);
slowFilt = filter(Hd, slow);

close all;
figure;
hold on;
if(false)
plot(x, fast);
plot(x, fastFilt);
else
plot(x, slow);
plot(x, slowFilt);
end
