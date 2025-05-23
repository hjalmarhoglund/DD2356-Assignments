#!/bin/bash

for i in {0..9};
do
	j=$((i*10));
	h1=$(sha256sum "wave_output_$j.txt" | awk '{print $1}');
	h2=$(sha256sum "ser_wave_output_$j.txt" | awk '{print $1}');
	h3=$(sha256sum "halo_wave_output_$j.txt" | awk '{print $1}');
	if [[ "$h1" == "$h2" && "$h2" == "$h3" ]];
	then
		echo "$j check OK!";
	else
		echo "$j failed! $h1 : $h2 : $h3";
	fi
done
