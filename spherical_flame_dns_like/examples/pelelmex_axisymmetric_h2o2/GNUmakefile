# PeleLMeX case template for a 2D RZ premixed H2/O2 spherical hot-kernel case.
# Set these paths here or export them in the shell before running make.
PELE_HOME ?= $(HOME)/PeleLMeX

DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE
USE_SYCL        = FALSE

Eos_Model       := Fuego
Transport_Model := Simple
Chemistry_Model := LiDryer

include $(PELE_HOME)/Exec/Make.PeleLMeX
