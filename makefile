NAME = KS.c

prova: $(NAME) 
	gcc $(NAME) -O3 -g -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -o prog