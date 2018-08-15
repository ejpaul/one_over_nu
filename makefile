LIBSTELL_DIR = mini_libstell
LIBSTELL = mini_libstell/mini_libstell.a

ifdef NERSC_HOST
        HOSTNAME = $(NERSC_HOST)
endif

ifeq ($(HOSTNAME),edison)
	FC = ftn
else ifeq ($(HOSTNAME),cori)
	FC = ftn
	EXTRA_COMPILE_FLAGS = -qopenmp 
	EXTRA_LINK_FLAGS =  -qopenmp 
else
	FC = mpif90
	EXTRA_COMPILE_FLAGS = -fopenmp -I/opt/local/include -ffree-line-length-none
    EXTRA_LINK_FLAGS =  -fopenmp -L/opt/local/lib -lnetcdff  -lnetcdf
endif

TARGET = one_over_nu

export

.PHONY: all clean

all: $(TARGET)

include makefile.depend

%.o: %.f90 $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

%.o: %.f $(LIBSTELL_DIR)/mini_libstell.a
	$(FC) $(EXTRA_COMPILE_FLAGS) -I $(LIBSTELL_DIR) -c $<

lib$(TARGET).a: $(OBJ_FILES)
	ar rcs lib$(TARGET).a $(OBJ_FILES)

$(TARGET): lib$(TARGET).a $(TARGET).o $(LIBSTELL)
	$(FC) -o $(TARGET) $(TARGET).o lib$(TARGET).a $(LIBSTELL) $(EXTRA_LINK_FLAGS)

$(LIBSTELL_DIR)/mini_libstell.a:
	$(MAKE) -C mini_libstell

clean:
	rm -f *.o *.mod *.MOD *~ $(TARGET) *.a
	cd mini_libstell; rm -f *.o *.mod *.MOD *.a
