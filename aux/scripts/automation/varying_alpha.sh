for alpha in $(seq 0.5 0.01 1.01); do
    echo "alpha = $alpha"
    ./build/main $alpha
    cp ./log/all-kmers.log ./results/alpha_exploration/$alpha.log
done