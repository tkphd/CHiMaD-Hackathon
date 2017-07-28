#!/bin/bash

if [[ -f ../square/electrochem.000100.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.000100.dat square.005.png
fi

if [[ -f ../square/electrochem.000200.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.000200.dat square.010.png
fi

if [[ -f ../square/electrochem.000400.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.000400.dat square.020.png
fi

if [[ -f ../square/electrochem.001000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.001000.dat square.050.png
fi

if [[ -f ../square/electrochem.002000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.002000.dat square.100.png
fi

if [[ -f ../square/electrochem.004000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.004000.dat square.200.png
fi

if [[ -f ../square/electrochem.008000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.008000.dat square.400.png
fi

if [[ -f ../square/electrochem.020000.dat ]]
then
	mmsp2png --field=0 ../square/electrochem.020000.dat square.1000.png
fi
