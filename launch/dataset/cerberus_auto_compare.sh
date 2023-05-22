#!/bin/bash
# for a given list of datasets and durations run vilo_autotest.sh 

ws_root='/home/shuoyang/Documents/vilo_dev/vilo_ws'
bag_folder='230517_sqh_mocap'

# a string of dataset names: 20230517_risqh_02speed_mocap, 20230517_risqh_04speed_mocap, ..
declare -a dataset_name_list=("20230517_risqh_02speed_mocap" 
                              "20230517_risqh_04speed_mocap"
                              "20230517_risqh_06speed_mocap" 
                              "20230517_risqh_08speed_mocap_more_turns" 
                              "20230517_risqh_10speed_mocap" )
declare -a dataset_duration_list=("80" 
                                  "44"
                                  "33" 
                                  "48" 
                                  "29")

# Iterate the string array using for loop           

for i in "${!dataset_name_list[@]}"; do
  printf ' %s, %s\n'  "${dataset_name_list[i]}"  "${dataset_duration_list[i]}"

    dataset_name=${dataset_name_list[i]}
    dataset_duration=${dataset_duration_list[i]}

    bag_path="$ws_root/bags/$bag_folder/$dataset_name.bag"
    ls $bag_path

    # run comparison
    temp_config="/tmp/temp_run_cerberus.yaml"
    cp $ws_root/src/Cerberus/config/go1_config/go1_realsense_left.yaml /tmp/go1_realsense_left.yaml
    cp $ws_root/src/Cerberus/config/go1_config/go1_realsense_right.yaml /tmp/go1_realsense_right.yaml
    cp $ws_root/src/Cerberus/config/go1_config/hardware_go1_vilo_config.yaml "$temp_config"

    sed -i "s/\(dataset_name: \)\"lab\"/\1\"$dataset_name\"/" "$temp_config"
    roslaunch vilo cerberus_auto.launch dataset_name:="$dataset_name" dataset_duration:="$dataset_duration"

sed -i 's/\(optimize_leg_bias: \)0/\11/' "$temp_config"
    roslaunch vilo cerberus_auto.launch dataset_name:="$dataset_name" dataset_duration:="$dataset_duration"

    # finally remove the temp file
    rm "$temp_config"
done