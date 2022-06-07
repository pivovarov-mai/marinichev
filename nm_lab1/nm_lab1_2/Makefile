CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wextra -Wall -Wno-sign-compare -Wno-unused-result
FILES = main.cpp
NAME = solution

all: solution

solution:
	$(CXX) $(CXXFLAGS) -o $(NAME) $(FILES)

debug:
	$(CXX) $(CXXFLAGS) -o $(NAME) $(FILES) -g3

valgrind:
	valgrind --leak-check=full --show-leak-kinds=all --log-file=valgrind.txt ./$(NAME) ./tests/test_02.txt

run:
	./$(NAME) ./tests/test_02.txt

clean:
	rm -f *.txt *.o $(NAME)