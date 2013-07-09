FLAGS = -Wall -std=c99

all:
	make 1d
	make 2d
	make level-n-trees

1d: 1d.c
	gcc $(FLAGS) -o ./1d.exe ./1d.c

2d: 2d.c
	gcc $(FLAGS) -o ./2d.exe ./2d.c

mt: mt.c
	gcc $(FLAGS) -o ./mt.o ./mt.c

level-n-trees: level-n-trees.c mt.o
	gcc $(FLAGS) -o ./level-n-trees.exe ./mt.o ./level-n-trees.c

clean:
	rm -f ./*.exe
	rm -f ./*.o
