# - Try to find rdkit
# Once done this will define
#
#  rdkit_FOUND          - True if rdkit found.
#  rdkit_INCLUDE_DIRS   - where to find the rdkit include directories
#  rdkit_LIBRARIES      - List of libraries when using rdkit.
#  rdkit_VERSION_STRING - The version of rdkit found (x.y.z)

INCLUDE( LibFindMacros)

# this package
SET( rdkit_PACKAGE_NAME rdkit)
#SET( rdkit_PACKAGE_VERSION 1.2022.09)
SET( rdkit_PACKAGE_VERSION ${rdkit_FIND_VERSION})

# Include directory
FIND_PATH(${rdkit_PACKAGE_NAME}_INCLUDE_SAMPLE NAMES GraphMol/Atom.h)
SET(${rdkit_PACKAGE_NAME}_INCLUDE_DIR ${${rdkit_PACKAGE_NAME}_INCLUDE_SAMPLE})
MESSAGE(STATUS "${rdkit_PACKAGE_NAME} top-level include directory: ${${rdkit_PACKAGE_NAME}_INCLUDE_DIR}")

# Names of individual libraries
SET( ${rdkit_PACKAGE_NAME}_LIBRARIES_ABREVS RDKitAbbreviations RDKitAlignment RDKitCatalogs RDKitChemicalFeatures RDKitChemReactions RDKitChemTransforms RDKitCIPLabeler RDKitcoordgen RDKitDataStructs RDKitDepictor RDKitDeprotect RDKitDescriptors RDKitDistGeometry RDKitDistGeomHelpers RDKitEigenSolvers RDKitFileParsers RDKitFilterCatalog RDKitFingerprints RDKitFMCS RDKitForceFieldHelpers RDKitForceField RDKitFragCatalog RDKitga RDKitGenericGroups RDKitGraphMol RDKithc RDKitInfoTheory RDKitmaeparser RDKitMMPA RDKitMolAlign RDKitMolCatalog RDKitMolChemicalFeatures RDKitMolDraw2D RDKitMolEnumerator RDKitMolHash RDKitMolInterchange RDKitMolStandardize RDKitMolTransforms RDKitO3AAlign RDKitOptimizer RDKitPartialCharges RDKitRDBoost RDKitRDGeneral RDKitRDGeometryLib RDKitRDStreams RDKitReducedGraphs RDKitRGroupDecomposition RDKitRingDecomposerLib RDKitScaffoldNetwork RDKitShapeHelpers RDKitSimDivPickers RDKitSLNParse RDKitSmilesParse RDKitSubgraphs RDKitSubstructLibrary RDKitSubstructMatch RDKitTautomerQuery RDKitTrajectory)

# Libraries
#MESSAGE(STATUS "All library abbreviations: ${${rdkit_PACKAGE_NAME}_LIBRARIES_ABREVS}")
FOREACH(X IN LISTS ${rdkit_PACKAGE_NAME}_LIBRARIES_ABREVS)
  FIND_LIBRARY( ${rdkit_PACKAGE_NAME}_${X}_LIBRARY NAMES ${X})
  #MESSAGE(STATUS "${rdkit_PACKAGE_NAME} current library: ${${rdkit_PACKAGE_NAME}_${X}_LIBRARY}")
  LIST(APPEND ${rdkit_PACKAGE_NAME}_LIBRARY_LIST ${${rdkit_PACKAGE_NAME}_${X}_LIBRARY})
ENDFOREACH()
SET(${rdkit_PACKAGE_NAME}_LIBRARIES ${${rdkit_PACKAGE_NAME}_LIBRARY_LIST})
MESSAGE(STATUS "List of ${rdkit_PACKAGE_NAME} libraries: ${${rdkit_PACKAGE_NAME}_LIBRARIES}")

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
#SET( ${rdkit_PACKAGE_NAME}_PROCESS_INCLUDES ${rdkit_PACKAGE_NAME}_INCLUDE_DIRS)
SET( ${rdkit_PACKAGE_NAME}_PROCESS_INCLUDES ${rdkit_PACKAGE_NAME}_INCLUDE_DIR)
SET( ${rdkit_PACKAGE_NAME}_PROCESS_LIBS ${rdkit_PACKAGE_NAME}_LIBRARIES)
LIBFIND_PROCESS( ${rdkit_PACKAGE_NAME})

# runtime library
#IF( MINGW)
#  FIND_LIBRARY( ${rdkit_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES rdkit1.dll)
#ELSEIF( APPLE)
#  FIND_LIBRARY( ${rdkit_PACKAGE_NAME}_RUNTIME_LIBRARY NAMES rdkit.dylib)
#ELSEIF( UNIX)
#IF(UNIX)
#    INCLUDE( MacroGetSoName)
#    SET( ${rdkit_PACKAGE_NAME}_RUNTIME_LIBRARY ${${rdkit_PACKAGE_NAME}_LIBRARIES} CACHE INTERNAL "runtime library for ${rdkit_PACKAGE_NAME}")
#    MACRO_GET_SO_NAME( ${${rdkit_PACKAGE_NAME}_LIBRARIES})
#    SET( ${rdkit_PACKAGE_NAME}_RUNTIME_LIBRARY_RENAME ${SO_FILE_NAME} CACHE INTERNAL "so name of runtime library for ${rdkit_PACKAGE_NAME}")
#    UNSET( SO_FILE_NAME)
#ENDIF()
