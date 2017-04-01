main:main.o
	gcc $< -o $@
main.o:main.c newCG.c
	gcc -c $< -o $@
clean: 
	rm *.o  main
