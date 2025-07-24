cmake_minimum_required(VERSION 3.12)
project(exercises)

file(GLOB_RECURSE EXERCISE_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach(EX_FILE ${EXERCISE_SOURCES})
    file(RELATIVE_PATH REL_PATH ${CMAKE_CURRENT_SOURCE_DIR} ${EX_FILE})

    # Sanitize for target name
    string(REPLACE "/" "_" EX_TARGET ${REL_PATH})
    string(REPLACE " " "_" EX_TARGET ${EX_TARGET})
    string(REPLACE "." "_" EX_TARGET ${EX_TARGET})
    string(REPLACE ".cpp" "" EX_TARGET ${EX_TARGET})

    add_executable(${EX_TARGET} ${EX_FILE})

    set_target_properties(${EX_TARGET} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/exercises"
    )
endforeach()
