set(util_source_files
  dictionary.h
  IO/ASCIItoBinary
  IO/ASCIItoBinary/LoadBinary.h
  IO/ASCIItoBinary/SaveBinary.h
  IO/Load.h
  IO/MPI_Save.h
  IO/MPI_TecplotOut.h
  IO/NewLoad.h
  MPI_Breakdown.h
  parse.h
  PrecisionTimer.h
  PrecisionTimerLinux.h
  Singleton.h
  Solvers/GMRES.h
  Solvers/GMRESbase.h
  Solvers/MPI_GMRES.h
  split.h
  Tensor.h
  TensorMath.h
  TensorMPI.h
  TensorStream.h
  Timer.h
  tostring.h
  VectorStreamPrint.h
)
add_library(util ${util_source_files})

SET_TARGET_PROPERTIES(util PROPERTIES LINKER_LANGUAGE C)
