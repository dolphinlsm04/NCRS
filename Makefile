LIBS = -D_DEFAULT_SOURCE -D_XOPEN_SOURCE=600 -lncursesw -ltinfo -ldrm -lm

ncrs : main.o inputbox.o
	gcc -o ncrs main.o inputbox.o $(LIBS)

inputbox.o : inputbox_ncrs.c
	gcc -c -o inputbox.o inputbox_ncrs.c
main.o : main.c
	gcc -c -o main.o main.c
	
clean:
	rm *.o ncrs
