cp -a bibli_fonctions $HOME
Ensuite créer un fichier nommé ccc contenant les deux ligne suivantes :
g++ -lm -Wall -I$HOME/bibli_fonctions $1 $HOME/bibli_fonctions/bibli_fonctions.ar && ./a.out
rm ./a.out 
et placer ce fichier dans le répertoire bin

