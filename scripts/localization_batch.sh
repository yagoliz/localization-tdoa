#!/usr/bin/zsh
#sudo rmmod dvb_usb_rtl28xxu rtl2832

# Iterations
num_iterations=58
start_num=142

if [ "$#" -ne 7 ]; then
  echo "Usage: $0 <directory-of-rtlsdr> <first-frequency (MHz)> <second-frequency (MHz)> <num-samples> <sampling-frequency> <device-num> <output-dir>"
  exit 1
fi

directory="$7"
current_day=$(date +%Y%m%d)

if [ ! -d "$directory/ltess/$2_$3_$4" ]
then
    mkdir -p "$directory/ltess/$2_$3_$4"
    echo "Directory: $directory/ltess/$2_$3_$4"
fi

if [ ! -d "$directory/localization/$2_$3_$4" ]
then
    mkdir -p "$directory/localization/$2_$3_$4"
fi

directory_ltess="$directory/ltess/$2_$3_$4"
directory_data="$directory/localization/$2_$3_$4"

# Create array for experiments
declare -a dates
for i in {1..$num_iterations}
do
    next_exp_time=$((i*15))
    dates[i]=$(date -d "$next_exp_time seconds" +%s.%N)
done

for i in {1..$num_iterations}
do
    echo "Loop: $i"
    # We get the time when we launch the librtlsdr 2 freqency
    echo ${dates[i]}

    # Launch regular RTL-SDR to get LTE samples for LO corretion
    rtl_sdr -f 806e6 -n 3840000 -s 1920000 -g 20 -d $6 $directory_ltess/D$6-E$((i+start_num))-ltess.dat

    # Testing purposes
    date_now=$(date +%s.%N)
    echo $date_now

    # Get time until when to sleep
    sleep_time=$((dates[i]-date_now))
    echo $sleep_time
    sleep $sleep_time

    echo "Command: $1/rtl_sdr -f $2e6 -h $3e6 -n $4 -s $5 -d $6 -g 20 $directory_data/D$6-E$((i+start_num))-localization.dat"
    echo "$(date +%s.%N)" | tee -a device_$6.txt
    $1/rtl_sdr -f $2e6 -h $3e6 -n $4 -s $5 -d $6 -g 20 $directory_data/D$6-E$((i+start_num))-localization.dat
done
