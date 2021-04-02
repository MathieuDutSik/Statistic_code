PROGRAM		= CorrEntropyFit

OBJECTS		= $(PROGRAM:%=%.o)

CC		= g++ -std=c++17 -Wall -O3 -g
#CC		= g++ -std=c++17 -Wall -O3 -g -fsanitize=address -fsanitize=pointer-subtract -fsanitize=undefined -fsanitize=null -fsanitize=bounds -fsanitize=bounds-strict -fsanitize=alignment

NOLINKFLAG 	= -c
OUTPUTFLAG 	= -o

CLOCAL          = -I../src_basic -I../basic_common_cpp   -I../src_number -I../src_matrix -I../src_comb -I../src_graph -I../src_latt -I../src_poly -I../src_perfect -I.


CFLAGS 		= -I$(GMP_INCDIR) -I$(PERMLIB_INCDIR) $(PERMUTALIB_INCLUDE) -I$(LIBBLISS_INCDIR) -I$(BOOST_INCDIR) -I$(EIGEN_PATH) $(MPI_INCLUDE_CPP) $(CDDLIB_INCLUDE) $(GLPK_INCLUDE) $(NAUTY_INCLUDE) $(NETCDF_CXX_ALLINC)


LIBS		= -lm $(GMP_CXX_LINK) -pthread $(LIBBLISS_LINK) $(BOOST_LINK) $(MPI_LINK_CPP) $(CDDLIB_LINK) $(GLPK_LINK) $(NAUTY_LINK) $(NETCDF_CXX_LINK)


%.o: %.cpp
	@echo
	@echo compile $<
	$(CC) $(OUTPUTFLAG) $@ $(CLOCAL) $(CFLAGS) $(NOLINKFLAG) $<

all:	$(PROGRAM)

$(PROGRAM): 	$(OBJECTS)
		@echo
		@echo link $@
		$(CC) $(LINKOBJECTS) $(DEBUGFLAG) $@.o $(LIBS) -o $@

clean:;
	rm -f *.o
	rm -f $(PROGRAM)
