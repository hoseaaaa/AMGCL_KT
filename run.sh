#!/bin/bash

#  运行格式
#  ./run.sh $1 $2 $3 
# $1 代表不同需求，amgcl 或具体项目
# $2 表示数据
# $3 是否需要make

set -e
time=$(date "+%Y%m%d-%H%M%S")

HOME_PATH=/home/yhz/AMGCL_KT
DATA_PATH=$HOME_PATH/amgcl_data
INFO_OUT=$HOME_PATH/info_out
POISSON3DB=$DATA_PATH/poisson3Db

TARGET_PATH=$HOME_PATH/build

data_G=$DATA_PATH/$2/ibmpg.G.txt
data_U=$DATA_PATH/$2/ibmpg.u.txt
data_outx=$DATA_PATH/$2/data_outx.txt



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
    "dc_keti")
        echo "Running $1... $2 // $3 "
        cd "$TARGET_PATH/tutorial/6.power_grid_dc"
        case "$3" in
            "make")      
                make -j8
            ;;
        esac
        echo "$time" > "${INFO_OUT}/6.power_grid_dc.txt"
        ./dc_keti  $data_G  $data_U $data_outx >> "${INFO_OUT}/6.power_grid_dc.txt"
    ;;
    "amgcl")
        echo "Running AMGCL..."
        rm -rf build/
        cmake -Bbuild -DAMGCL_BUILD_TESTS=ON -DAMGCL_BUILD_EXAMPLES=ON .
    ;;
    "all")
        echo "Running $1..all ."
        cd "$TARGET_PATH/tutorial/6.power_grid_dc"
        case "$3" in
            "make")      
                make -j8
            ;;
        esac
        rm -rf  $HOME_PATH/info_out/*
        echo "$time" > "${INFO_OUT}/6.power_grid_dc.txt"
        ./dc_keti  $HOME_PATH/amgcl_data/ibmpg1/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg1/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg1/data_outx.txt >> "${INFO_OUT}/01.ibmpg.info" || true
        ./dc_keti  $HOME_PATH/amgcl_data/ibmpg2/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg2/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg2/data_outx.txt >> "${INFO_OUT}/02.ibmpg.info" || true
        ./dc_keti  $HOME_PATH/amgcl_data/ibmpg3/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg3/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg3/data_outx.txt >> "${INFO_OUT}/03.ibmpg.info" || true
        ./dc_keti  $HOME_PATH/amgcl_data/ibmpg4/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg4/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg4/data_outx.txt >> "${INFO_OUT}/04.ibmpg.info" || true
        ./dc_keti  $HOME_PATH/amgcl_data/ibmpg5/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg5/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg5/data_outx.txt >> "${INFO_OUT}/05.ibmpg.info" || true
        ./dc_keti  $HOME_PATH/amgcl_data/ibmpg6/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg6/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg6/data_outx.txt >> "${INFO_OUT}/06.ibmpg.info" || true
    ;;
esac