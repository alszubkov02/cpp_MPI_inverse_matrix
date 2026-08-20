CFLAGS = -isystem /usr/lib/x86_64-linux-gnu/openmpi/include -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format

a.out:	main.o func.o solve.o
		mpicxx $^ $(CFLAGS) -o $@
main.o: main.cpp func.h
		mpicxx -c $(CFLAGS) main.cpp
func.o: func.cpp func.h
		mpicxx -c $(CFLAGS) func.cpp
solve.o: solve.cpp func.h
		mpicxx -c $(CFLAGS) solve.cpp	
clean:
		rm -f main.o func.o solve.o a.out
