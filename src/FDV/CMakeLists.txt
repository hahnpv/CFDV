set(FDV_source_files
  FDVGalerkin.h
  Flow.h
  NavierStokes2D.h
  Node.h
  Element.h
  NavierStokes3D.h
  tau.h
  Flow_FGB.h
  NavierStokes1D.h
  NavierStokes.h
  Thermo.h
  Function/ApplyBC.h
  Function/ClearElement.h
  Function/eval2dboundary.h
  Function/MatrixOut.h
  Function/NodeCheck.h
  Function/AxisymmetricFlow.h
  Function/EnforceBC.h
  Function/FDVParam.h
  Function/MPI_CFL.h
  Function/NodeUnpack.h
  Function/CalcLength.h
  Function/eval1dboundary.h
  Function/LUFactorization.h
  Function/MPI_RMSError.h
  Function/RMSBase.h
)
add_library(FDV ${FDV_source_files})

SET_TARGET_PROPERTIES(FDV PROPERTIES LINKER_LANGUAGE C)
