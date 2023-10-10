#!/bin/bash
set -e
time=$(date "+%Y%m%d-%H%M%S")

HOME_PATH=/home/yhz/amgcl
DATA_PATH=$HOME_PATH/amgcl_data
INFO_OUT=$HOME_PATH/info_out
POISSON3DB=$DATA_PATH/poisson3Db

TARGET_PATH=$HOME_PATH/build


data_C=$DATA_PATH/ibmpg3t_12456_no_Vgnd/CSC_C.txt
data_G=$DATA_PATH/ibmpg3t_12456_no_Vgnd/CSC_G.txt
data_B=$DATA_PATH/ibmpg3t_12456_no_Vgnd/CSC_B.txt
data_U=$DATA_PATH/ibmpg3t_12456_no_Vgnd/u_t.txt
data_outx=$DATA_PATH/ibmpg3t_12456_no_Vgnd/data_outx.txt

case "$1" in
    "poisson3Db")
        echo "Running poisson3Db..."
        cd "$TARGET_PATH/tutorial/1.poisson3Db"
        case "$2" in
            "make")      
                make -j8
            ;;
        esac
        echo "$time" > "${INFO_OUT}/1.POISSON3DBinfo.txt"
        ./poisson3Db "${POISSON3DB}/poisson3Db.mtx" "${POISSON3DB}/poisson3Db_b.mtx" >> "${INFO_OUT}/1.POISSON3DBinfo.txt"
    ;;
    "dc_test")
        echo "Running dc_test..."
        cd "$TARGET_PATH/tutorial/6.power_grid_dc"
        case "$2" in
            "make")      
                make -j8
            ;;
        esac
        echo "$time" > "${INFO_OUT}/6.power_grid_dc.txt"
        ./dc_test $data_C $data_G $data_B $data_U $data_outx >> "${INFO_OUT}/6.power_grid_dc.txt"
    ;;
    "amgcl")
        echo "Running AMGCL..."
        rm -rf build/
        cmake -Bbuild -DAMGCL_BUILD_TESTS=ON -DAMGCL_BUILD_EXAMPLES=ON .
    ;;

esac