default:
	gcc -lomp -Xclang -fopenmp -O3 main.c
	./a.out
createplot:
	cat out_c1.txt out_c32.txt out_c64.txt out_c128.txt | python3 plot_data.py
