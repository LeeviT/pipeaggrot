cmake_minimum_required(VERSION 3.14)
project(cppsrc)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(cppsrc
        Combo_class.cpp
        Combo_class.hpp
        dist_combo.hpp
        main.cpp
        main_velo.cpp
        makefile
        orientation_control_vol.hpp
        orientation_control_vol_acorrected.cpp
        runscript_layout.sh
        triple.hpp)
