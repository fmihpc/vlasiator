# Standalone: needs only a C++20 compiler, no VLASIATOR_ARCH.
CXX ?= g++
FLAGS = -W -Wall -Wextra -pedantic -std=c++20 -O2 -I../..

default: test_erase

test_erase: test_erase.cpp ../../open_bucket_hashtable.h
	${CXX} ${FLAGS} $< -o $@

check: test_erase
	./test_erase

clean:
	rm -f test_erase
