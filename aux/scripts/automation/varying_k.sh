for K in $(seq 15 1 32); do
    echo "k = $K"
    ./build/main 0.5 $K
    python3 nothing.py -k $K
done