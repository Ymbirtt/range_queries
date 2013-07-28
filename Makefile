FLAGS = -Wall -g -O3 --std=c99 -D__USE_MINGW_ANSI_STDIO
CC = gcc

all:
	make 1d
	make 2d
	make level-n-trees
	make nd-trees

1d: 1d.c
	$(CC) $(FLAGS) -o ./1d.exe ./1d.c

2d: 2d.c
	$(CC) $(FLAGS) -o ./2d.exe ./2d.c

mt: mt.c
	$(CC) $(FLAGS) -o ./mt.o ./mt.c

level-n-trees: level-n-trees.c mt.o
	$(CC) $(FLAGS) -o ./level-n-trees.exe ./mt.o ./level-n-trees.c

nd-trees: nd-trees.c mt.o
	$(CC) $(FLAGS) -o ./nd-trees.exe ./mt.o ./nd-trees.c

clean:
	rm -f ./*.exe
	rm -f ./*.o
