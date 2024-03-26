CC=-c -Wall --pedantic
bibli_perso.ar : affichage.o matrices.o 
	ar -r bibli_perso.ar affichage.o matrices.o 
affichage.o: affichage.cpp
	g++ -g -I/public/mphyo/bibli_fonctions -I/usr/include/python3.10 /public/mphyo/bibli_fonctions/bibli_fonctions.ar -lm -lpython3.10 $(CC) affichage.cpp
matrices.o : matrices.cpp
	g++ -g -I/public/mphyo/bibli_fonctions -I/usr/include/python3.10 /public/mphyo/bibli_fonctions/bibli_fonctions.ar -lm -lpython3.10 $(CC) matrices.cpp