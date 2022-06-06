CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wextra -Wall -Wno-sign-compare -Wno-unused-result
FILES = main.cpp
NAME = generator
ROWS = 100
COLS = 100

all: generator

generator:
	$(CXX) $(CXXFLAGS) -o $(NAME) $(FILES)

debug:
	$(CXX) $(CXXFLAGS) -o $(NAME) $(FILES) -g3

valgrind:
	valgrind --leak-check=full --show-leak-kinds=all --log-file=valgrind.txt ./$(NAME) $(ROWS) $(COLS)

run:
	./$(NAME) $(ROWS) $(COLS)

clean:
	rm -f *.txt *.o $(NAME)