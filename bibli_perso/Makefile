CC=-c -Wall --pedantic
bibli_perso.ar : affichage.o matrices.o 
	ar -r bibli_perso.ar affichage.o matrices.o 
affichage.o: affichage.cpp
	g++ $(CC) affichage.cpp
matrices.o : matrices.cpp
	g++ $(CC) matrices.cpp