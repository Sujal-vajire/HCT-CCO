# HCT-CCO build configuration (modified from OpenCCO). See top-level README and LICENSE.md.
# Modifications Copyright (C) 2024-2026 S. L. Vajire, J. S. Choy, G. S. Kassab and L.-C. Lee;
# original Copyright (C) 2023 B. Kerautret, Phuc Ngo, N. Passat, H. Talbot, C. Jaquet. GNU GPL v3.

SET(BIN_SRC
  HctCco
)

FOREACH(FILE ${BIN_SRC})
  add_executable(${FILE} ${FILE}.cpp ${SOURCE_BASE_FILES} ${SOURCE_BASE_HEADER_FILES} )
  target_link_libraries (${FILE} ${CERES_LIBRARIES} ${DGTAL_LIBRARIES} )
  # Path to the consolidated post-processing script, invoked automatically after
  # generation (disable with --no-postprocess).
  target_compile_definitions(${FILE} PRIVATE
    POSTPROCESS_SCRIPT="${CMAKE_SOURCE_DIR}/postprocessing/coronary_postprocess.py")
ENDFOREACH(FILE)
