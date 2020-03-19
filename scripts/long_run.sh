#!/usr/bin/zsh
#sudo rmmod dvb_usb_rtl28xxu rtl2832


if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <device-num> <output-dir>"
  exit 1
fi

directory="$2"
current_day=$(date +%Y%m%d)

if [ ! -d "$directory/ppm_stability" ]
then
    mkdir -p "$directory/ppm_stability"
    echo "Directory: $directory/ppm_stability"
fi

# Create directory
directory_ppm="$directory/ppm_stability"

# Sleep until 10 minutes and 15 seconds
test_time=$(date -d "10 minutes 15 seconds" +%s.%N)
echo $test_time

# First we heat up the RTL-SDR
timeout 10m rtl_test -d $1

# Launch regular RTL-SDR to get LTE samples for LO corretion
rtl_sdr -f 806e6 -n 3840000 -s 1920000 -g 20 -d $1 $directory_ppm/D$1_ltess.dat

# Get time until when to sleep
date_now=$(date +%s.%N)
sleep_time=$((test_time-date_now))
echo $sleep_time
sleep $sleep_time

rtl_sdr -f 806e6 -n 200e6 -s 2e6 -g 20 -d $1 $directory_ppm/D$1_ppm.dat
