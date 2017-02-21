function out = coinflip0(numFlips)

heads = 0;
for i = 1:numFlips
    r = rand;
    if rand <= 0.5
        heads = heads + 1;
    end
end
out = heads/numFlips;
