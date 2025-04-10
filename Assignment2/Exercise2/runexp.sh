#!/bin/bash

for i in {0..9}; do
    ./a1e6.out > "res1e6_$i.txt";
done

for i in {0..9}; do
    ./a1e7.out > "res1e7_$i.txt";
done

for i in {0..9}; do
    ./a1e8.out > "res1e8_$i.txt";
done

for i in {0..9}; do
    ./a1e9.out > "res1e9_$i.txt";
done
