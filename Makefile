FLAGS = -Wall -std=c99

all:
	make 1d
	make 2d
	make nd

1d: 1d.c
	gcc $(FLAGS) -o ./1d.exe ./1d.c

2d: 2d.c
	gcc $(FLAGS) -o ./2d.exe ./2d.c

nd: nd.c
	gcc $(FLAGS) -o ./nd.exe ./nd.c

clean:
	rm -rf ./*.exe
