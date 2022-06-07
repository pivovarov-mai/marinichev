CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wextra -Wall -Wno-sign-compare -Wno-unused-result
FILES = main.cpp
NAME = solution
PRECISION = 0.001
METHOD = -n

all: solution

solution:
	$(CXX) $(CXXFLAGS) -o $(NAME) $(FILES)

debug:
	$(CXX) $(CXXFLAGS) -o $(NAME) $(FILES) -g3

valgrind:
	valgrind --leak-check=full --show-leak-kinds=all --log-file=valgrind.txt ./$(NAME) ./tests/100x100.txt -p $(PRECISION) $(METHOD)

run:
	./$(NAME) ./tests/100x100.txt -p $(PRECISION) $(METHOD)

clean:
	rm -f *.txt *.o $(NAME)