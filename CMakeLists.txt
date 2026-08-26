cmake_minimum_required(VERSION 3.16)

project(Warper VERSION 1.0.0 LANGUAGES CXX)

set(SUPERCOLLIDER_SOURCE_DIR "" CACHE PATH
    "Path to the SuperCollider source tree or SDK installation")

find_path(SC_PLUGIN_INTERFACE_DIR
    NAMES SC_PlugIn.h
    HINTS
        "${SUPERCOLLIDER_SOURCE_DIR}/include/plugin_interface"
        "${SUPERCOLLIDER_SOURCE_DIR}/include"
)

find_path(SC_COMMON_DIR
    NAMES SC_Types.h
    HINTS
        "${SUPERCOLLIDER_SOURCE_DIR}/include/common"
        "${SUPERCOLLIDER_SOURCE_DIR}/include"
)

if(NOT SC_PLUGIN_INTERFACE_DIR OR NOT SC_COMMON_DIR)
    message(FATAL_ERROR
        "SuperCollider SDK headers not found. Configure with "
        "-DSUPERCOLLIDER_SOURCE_DIR=/path/to/supercollider")
endif()

add_library(PolyFreqWarp MODULE PolyFreqWarp.cpp)

target_compile_features(PolyFreqWarp PRIVATE cxx_std_17)
target_include_directories(PolyFreqWarp PRIVATE
    "${SC_PLUGIN_INTERFACE_DIR}"
    "${SC_COMMON_DIR}"
)

set_target_properties(PolyFreqWarp PROPERTIES
    PREFIX ""
    OUTPUT_NAME "PolyFreqWarp"
    SUFFIX ".scx"
)

include(GNUInstallDirs)

if(APPLE AND NOT DEFINED SUPERCOLLIDER_EXTENSION_ROOT)
    set(SUPERCOLLIDER_EXTENSION_ROOT "$ENV{HOME}/Library/Application Support" CACHE PATH
        "Install path prefix, prepended onto install directories." FORCE)
endif()

set(SUPERCOLLIDER_EXTENSION_ROOT "SuperCollider/Extensions/Warper"
    CACHE PATH "Install destination for this SuperCollider extension package")
set(SUPERCOLLIDER_PLUGIN_DIR "${SUPERCOLLIDER_EXTENSION_ROOT}/plugins"
    CACHE PATH "Install destination for SuperCollider plugins")
set(SUPERCOLLIDER_EXTENSION_DIR "${SUPERCOLLIDER_EXTENSION_ROOT}/classes"
    CACHE PATH "Install destination for SuperCollider class files")
set(SUPERCOLLIDER_HELP_DIR "${SUPERCOLLIDER_EXTENSION_ROOT}/HelpSource/Classes"
    CACHE PATH "Install destination for SuperCollider help files")

install(TARGETS PolyFreqWarp
    LIBRARY DESTINATION "${SUPERCOLLIDER_PLUGIN_DIR}"
    RUNTIME DESTINATION "${SUPERCOLLIDER_PLUGIN_DIR}"
)

install(FILES PolyFreqWarp.sc
    DESTINATION "${SUPERCOLLIDER_EXTENSION_DIR}"
)

install(FILES HelpSource/Classes/PV_PolyFreqWarp.schelp
    DESTINATION "${SUPERCOLLIDER_HELP_DIR}"
)
