set(util_source_files
  IO/ASCIItoBinary
  IO/ASCIItoBinary/LoadBinary.h
  IO/ASCIItoBinary/SaveBinary.h
  IO/MPI_TecplotOut.h
  MPI_Breakdown.h
  parse.h
  PrecisionTimer.h
  Sensor.h
  Singleton.h
  Solvers/GMRES.h
  Solvers/GMRESbase.h
  Solvers/MPI_GMRES.h
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
