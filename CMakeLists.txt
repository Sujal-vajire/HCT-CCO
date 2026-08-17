# HCT-CCO build configuration, a coronary-specific extension of OpenCCO.
# Modifications Copyright (C) 2024-2026 S. L. Vajire, J. S. Choy, G. S. Kassab and
# L.-C. Lee (Michigan State University; California Medical Innovations Institute),
# modified from the original OpenCCO CMake files (Copyright (C) 2023 B. Kerautret,
# Phuc Ngo, N. Passat, H. Talbot, C. Jaquet). Released under the GNU GPL v3.
# See README for the list of changes.

cmake_minimum_required(VERSION 3.5)
project(HCT_CCO)
set (CMAKE_CXX_STANDARD 11)

# Header set shared by the generator (unchanged from the original implementation).
SET(SOURCE_BASE_HEADER_FILES
  ${PROJECT_SOURCE_DIR}/src/helpers/CoronaryGeometry.h
  ${PROJECT_SOURCE_DIR}/src/helpers/CoronaryGrowth.h
  ${PROJECT_SOURCE_DIR}/src/HctccoModel.h
  ${PROJECT_SOURCE_DIR}/src/helpers/CoronaryTreeIO.h
  ${PROJECT_SOURCE_DIR}/src/MyocardialDomain.h
)
SET(SOURCE_BASE_FILES
  ${PROJECT_SOURCE_DIR}/src/HctccoModel.ih
)

find_package(Ceres REQUIRED)

# CLI11 (bundled), local headers
include_directories( "${PROJECT_SOURCE_DIR}/ext/" "${PROJECT_SOURCE_DIR}/src/helpers/" "${PROJECT_SOURCE_DIR}/src/" )
FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS} )

add_subdirectory(bin)
