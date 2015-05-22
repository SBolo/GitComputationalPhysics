NAME = VMC.cpp

all: $(NAME) 
	g++ $(NAME) -O2 -lgsl -lgslcblas -larmadillo -I/opt/local/include -L/opt/local/lib -o prog

debug: $(NAME)
	g++ $(NAME) -O2 -g -lgsl -lgslcblas -larmadillo -I/opt/local/include -L/opt/local/lib -o prog
