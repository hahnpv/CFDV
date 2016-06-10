set(adap_source_files
  ElementAssociation.h
  MeshRefine2D.h
  MeshRefine.h
  MeshUnrefine2D.h
  RefineListAdaptive.h
  RefineListProgrammed.h
)
add_library(adap ${adap_source_files})
SET_TARGET_PROPERTIES(adap PROPERTIES LINKER_LANGUAGE C)
