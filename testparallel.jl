# Ordinary for loop.
a = 0;
@time for i = 1:20000
     a += Int(rand(Bool));
end

# Parallel for loop
@time @parallel (+) for i = 1:20000
    Int(rand(Bool));
end

# Predefined function.
@time sum(rand(0:1, 20000))

# For a not so small amount of work:
n = 200000000;
@time @parallel (+) for i = 1:n
    Int(rand(Bool))
end

# DO NOT run sum(rand(0:1, n)) or the ordinary for loop with 'n': You may run
# out of memory.
