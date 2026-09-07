FC = gfortran
FFLAGS = -O2 -fopenmp

SRCS = precision_mod.f90 linear_solver.f90 newton_module.f90 test_newton.f90
OBJS = $(SRCS:.f90=.o)
MODS = precision_mod.mod linear_solver.mod newton_module.mod
TARGET = newton_solver

all: $(TARGET)

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<
	
precision_mod.o: precision_mod.f90
linear_solver.o: linear_solver.f90 precision_mod.mod
newton_module.o: newton_module.f90 precision_mod.mod linear_solver.mod
test_newton.o: test_newton.f90 newton_module.mod

precision_mod.mod: precision_mod.o
linear_solver.mod: linear_solver.o
newton_module.mod: newton_module.o

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f *.o *.mod $(TARGET)
	
.PHONY: all clean
