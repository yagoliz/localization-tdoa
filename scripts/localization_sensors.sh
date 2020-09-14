#!/bin/bash
#sudo rmmod dvb_usb_rtl28xxu rtl2832

# Iterations
num_iterations=4
start_num=0

if [ "$#" -ne 5 ]; then
  echo "Usage: $0 <directory-of-rtlsdr> <first-frequency (MHz)> <second-frequency (MHz)> <num-samples> <sampling-frequency>"
  exit 1
fi

# Create array for experiments
declare -a dates
for i in $(seq 0 $num_iterations)
do
    next_exp_time=$(((i+1)*15))
    dates[i]=$(date -d "$next_exp_time seconds" +%s.%N)
    echo "${dates[i]}"
done

for i in $(seq 0 $num_iterations)
do
    echo "Loop: $i"

    # Launch regular RTL-SDR to get LTE samples for LO corretion
    rtl_sdr -f 806e6 -n 3840000 -s 1920000 -g 20 E$((i+start_num))-ltess.dat

    # Testing purposes
    date_now=$(date +%s.%N)

    # Get time until when to sleep
    sleep_time=`bc <<< "${dates[i]}-${date_now}"`
    echo $sleep_time
    sleep $sleep_time

    echo "Command: $1/rtl_sdr -f $2e6 -h $3e6 -n $4 -s $5 -d $6 -g 20 E$((i+start_num))-localization.dat"
    $1/rtl_sdr -f $2e6 -h $3e6 -n $4 -s $5 -g 30 E$((i+start_num))-localization.dat
done
