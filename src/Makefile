CFLAGS = -Wall -Werror -pedantic -std=c++17

default: curve_fit conj_grad precond_conj_grad

curve_fit: curve_fit.cpp
	g++ curve_fit.cpp $(CFLAGS) -o curve_fit

conj_grad: conj_grad.cpp
	g++ conj_grad.cpp $(CFLAGS) -o conj_grad

precond_conj_grad: precond_conj_grad.cpp
	g++ precond_conj_grad.cpp $(CFLAGS) -o precond_conj_grad

clean:
	-rm -f curve_fit
	-rm -f conj_grad
	-rm -f precond_conj_grad
	-rm -f "conj_grad data "*".txt"
	-rm -f "precond_conj_grad data "*".txt"
	-rm -f "data_"*".txt"