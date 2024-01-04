#!/bin/bash

#  运行格式
#  ./run.sh $1 $2 $3 
# $1 代表不同需求，amgcl 或具体项目
# $2 表示数据
# $3 是否需要make
set -e
time=$(date "+%Y%m%d-%H%M%S")

current_path=$(pwd)
export HOME_PATH="$current_path"
# HOME_PATH=/home/yhz/AMGCL_KT
DATA_PATH=$HOME_PATH/amgcl_data
INFO_OUT=$HOME_PATH/info_out
POISSON3DB=$DATA_PATH/poisson3Db

TARGET_PATH=$HOME_PATH/build

prefix=$(echo "$2" | cut -c1-5)
echo $prefix
data_G=$DATA_PATH/$2/$prefix.G.txt
data_U=$DATA_PATH/$2/$prefix.u.txt
data_outx=$DATA_PATH/$2/$prefix.AMG_x.txt



case "$1" in
    "make")
        echo "Running makefile make a..."
        cd "$TARGET_PATH/tutorial/6.power_grid_dc"
        make -j8 
    ;;
    "poisson3Db")
        echo "Running poisson3Db..."
        cd "$TARGET_PATH/tutorial/1.poisson3Db"
        echo "$time" > "${INFO_OUT}/1.POISSON3DBinfo.txt"
        ./poisson3Db "${POISSON3DB}/poisson3Db.mtx" "${POISSON3DB}/poisson3Db_b.mtx" >> "${INFO_OUT}/1.POISSON3DBinfo.txt"
    ;;
    "dc_keti")
        echo "Running $1... $2 // "
        cd "$TARGET_PATH/tutorial/6.power_grid_dc"
        echo "$time" > "${INFO_OUT}/00.$2.txt"
        ./dc_keti  $data_G  $data_U $data_outx >> "${INFO_OUT}/00.$2.txt"
    ;;
    "amgcl")
        echo "Running AMGCL..."
        rm -rf build/
        if [ -z "$2"]; then 
            cmake -Bbuild -DAMGCL_BUILD_TESTS=ON -DAMGCL_BUILD_EXAMPLES=ON .
        else
            cmake -Bbuild -DAMGCL_BUILD_TESTS=ON -DAMGCL_BUILD_EXAMPLES=ON  -DEPS_STRONG=$2 -DRELAX=$3 .
        fi
    ;;
    "all")
        echo "Running $1.$2.all ."
        cd "$TARGET_PATH/tutorial/6.power_grid_dc"
        case "$2" in 
            "ibm")
                echo "$time" > "${INFO_OUT}/6.$2.txt"
                ./dc_keti  $HOME_PATH/amgcl_data/ibmpg1/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg1/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg1/data_outx.txt > "${INFO_OUT}/i1.ibmpg.info" || true
                ./dc_keti  $HOME_PATH/amgcl_data/ibmpg2/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg2/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg2/data_outx.txt > "${INFO_OUT}/i2.ibmpg.info" || true
                ./dc_keti  $HOME_PATH/amgcl_data/ibmpg3/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg3/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg3/data_outx.txt > "${INFO_OUT}/i3.ibmpg.info" || true
                ./dc_keti  $HOME_PATH/amgcl_data/ibmpg4/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg4/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg4/data_outx.txt > "${INFO_OUT}/i4.ibmpg.info" || true
                ./dc_keti  $HOME_PATH/amgcl_data/ibmpg5/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg5/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg5/data_outx.txt > "${INFO_OUT}/i5.ibmpg.info" || true
                ./dc_keti  $HOME_PATH/amgcl_data/ibmpg6/ibmpg.G.txt  $HOME_PATH/amgcl_data/ibmpg6/ibmpg.u.txt $HOME_PATH/amgcl_data/ibmpg6/data_outx.txt > "${INFO_OUT}/i6.ibmpg.info" || true
            ;;
            "thu")
                echo "$time" > "${INFO_OUT}/6.$2.txt"
                ./dc_keti  $HOME_PATH/amgcl_data/thupg1/thupg.G.txt  $HOME_PATH/amgcl_data/thupg1/thupg.u.txt $HOME_PATH/amgcl_data/thupg1/data_outx.txt > "${INFO_OUT}/t1.thupg.info" || true
                ./dc_keti  $HOME_PATH/amgcl_data/thupg2/thupg.G.txt  $HOME_PATH/amgcl_data/thupg2/thupg.u.txt $HOME_PATH/amgcl_data/thupg2/data_outx.txt > "${INFO_OUT}/t2.thupg.info" || true
                ./dc_keti  $HOME_PATH/amgcl_data/thupg3/thupg.G.txt  $HOME_PATH/amgcl_data/thupg3/thupg.u.txt $HOME_PATH/amgcl_data/thupg3/data_outx.txt > "${INFO_OUT}/t3.thupg.info" || true
                ./dc_keti  $HOME_PATH/amgcl_data/thupg4/thupg.G.txt  $HOME_PATH/amgcl_data/thupg4/thupg.u.txt $HOME_PATH/amgcl_data/thupg4/data_outx.txt > "${INFO_OUT}/t4.thupg.info" || true
                # ./dc_keti  $HOME_PATH/amgcl_data/thupg5/thupg.G.txt  $HOME_PATH/amgcl_data/thupg5/thupg.u.txt $HOME_PATH/amgcl_data/thupg5/data_outx.txt > "${INFO_OUT}/t5.thupg.info" || true
                # ./dc_keti  $HOME_PATH/amgcl_data/thupg6/thupg.G.txt  $HOME_PATH/amgcl_data/thupg6/thupg.u.txt $HOME_PATH/amgcl_data/thupg6/data_outx.txt > "${INFO_OUT}/t6.thupg.info" || true
                # ./dc_keti  $HOME_PATH/amgcl_data/thupg7/thupg.G.txt  $HOME_PATH/amgcl_data/thupg7/thupg.u.txt $HOME_PATH/amgcl_data/thupg7/data_outx.txt > "${INFO_OUT}/t7.thupg.info" || true
                # ./dc_keti  $HOME_PATH/amgcl_data/thupg8/thupg.G.txt  $HOME_PATH/amgcl_data/thupg8/thupg.u.txt $HOME_PATH/amgcl_data/thupg8/data_outx.txt > "${INFO_OUT}/t8.thupg.info" || true
            ;;
        esac
    ;;
    "error")  
        cd "$HOME_PATH/error"
        ./error.sh
    ;;
    "clean")
        rm -rf  $HOME_PATH/info_out/*.info
    ;;
esac