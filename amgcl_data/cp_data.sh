current_path=$(pwd)
export CURR_PATH="$current_path"
IBM_MATRIX=$HOME/DCSIM_Improve/data_out/ibm_data_cq
THU_MATRIX=$HOME/DCSIM_Improve/data_out/thu_data_cq

cd $CURR_PATH  && rm -rf */
if [ "$1" = "ibm" ]; then
    cd "$IBM_MATRIX" && cp -r */ "$CURR_PATH"
elif [ "$1" = "thu" ]; then
    cd "$THU_MATRIX" && cp -r */ "$CURR_PATH"
elif [ "$1" = "all" ]; then
    cd "$IBM_MATRIX" && cp -r */ "$CURR_PATH"
    cd "$THU_MATRIX" && cp -r */ "$CURR_PATH"
else
    echo "未知选项: $1"
fi
