# Compilateur utilisé
CC=gfortran

# Le nom de l'exécutable
PROG = run

# Les fichiers source à compiler
SRC = Main.f90 Num.f90 Func.f90 Data.f90 System.f90 GC.f90

# La commande complète : compile seulement si un fichier a été modifié
$(PROG) : $(SRC)
	$(CC) $(SRC) -o $(PROG)
# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ $(PROG)
