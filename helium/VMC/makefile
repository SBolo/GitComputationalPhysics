NAME = vmc.cpp

all: $(NAME) 
	g++ $(NAME) -O2 -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -o prog

debug: $(NAME)
	g++ $(NAME) -g -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -o prog

thread: $(NAME)
	clang++ $(NAME) -lgsl -lgslcblas -I/opt/local/include -L/opt/local/lib -std=c++11 -stdlib=libc++ -o prog