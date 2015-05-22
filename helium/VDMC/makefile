NAME = vmc.cpp
NAME2 = dmc.cpp

vmc: $(NAME) 
	g++ $(NAME) -O2 -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -o vmc

dmc: $(NAME2)
	g++ $(NAME2) -O2 -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -o dmc

dmc_debug: $(NAME2)
	g++ $(NAME2) -g -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -o dmc

vmc_debug: $(NAME) 
	g++ $(NAME) -g -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -o vmc
