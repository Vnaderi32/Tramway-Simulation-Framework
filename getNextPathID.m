function id = getNextPathID(reset)
    persistent counter

    if nargin == 1 && reset
        counter = 0;
    end

    if isempty(counter)
        counter = 0;
    end

    id = counter;
    counter = counter + 1;
end