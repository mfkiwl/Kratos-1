set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosDemStructuresCouplingApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/DEMApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/PoromechanicsApplication)

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/DEMApplication )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication )

## generate variables with the sources
set(KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/dem_structures_coupling_application.cpp

    ## CONDITIONS
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/line_load_from_DEM_condition_2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/surface_load_from_DEM_condition_3d.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/multiaxial_control_module_fem_dem_generalized_2d_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/permeability_tensor_communicator_utility.cpp
)

set(KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/dem_structures_coupling_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/post_process_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/sand_production_utilities.cpp
)


add_library(KratosDemStructuresCouplingCore SHARED ${KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_CORE})
target_link_libraries(KratosDemStructuresCouplingCore PUBLIC KratosDEMCore KratosStructuralMechanicsCore KratosCore KratosPoromechanicsCore)
set_target_properties(KratosDemStructuresCouplingCore PROPERTIES COMPILE_DEFINITIONS "DEM_STRUCTURES_COUPLING_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosDemStructuresCouplingApplication MODULE THIN_LTO ${KRATOS_DEM_STRUCTURES_COUPLING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosDemStructuresCouplingApplication PUBLIC KratosDemStructuresCouplingCore KratosPoromechanicsCore)
set_target_properties(KratosDemStructuresCouplingApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosDemStructuresCouplingCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosDemStructuresCouplingApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

if(${ACTIVATE_DEBUG_MACRO} MATCHES ON) #MSI: Flag defined for debug Macro
    add_definitions(-DDEBUG_MACRO)
endif(${ACTIVATE_DEBUG_MACRO} MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosDemStructuresCouplingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosDemStructuresCouplingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/DemStructuresCouplingApplication.py KratosMultiphysics/DemStructuresCouplingApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosDemStructuresCouplingCore DESTINATION libs )
install(TARGETS KratosDemStructuresCouplingApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosDemStructuresCouplingCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosDemStructuresCouplingApplication" PARENT_SCOPE)