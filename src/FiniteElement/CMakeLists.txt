set(finite_element_source_files
  calc_TF_norm3D.h
  Face.h
  get_face.h
  TestFunctions1D.h
  TestFunctions3D.h
  TestFunctionsTriangle2D.h
  GaussQuad.h
  Lagrange.h
  TestFunctions2D.h
  TestFunctions.h
)
add_library(finite_element ${finite_element_source_files})

SET_TARGET_PROPERTIES(finite_element PROPERTIES LINKER_LANGUAGE C)
