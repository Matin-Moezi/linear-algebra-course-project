CC=gcc

main: main.c
	$(CC) main.c -o main -lm
deb: main.c
	$(CC) -g main.c -o deb -lm
