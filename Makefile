FLAGS = -Wall -std=c99

all: main.c
	gcc $(FLAGS) -o main.exe ./main.c

clean:
	rm -rf ./*.exe
