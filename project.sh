#!/usr/bin/bash

# # echo "File,Reciprocal,,Real" > parallel1.csv
# i=5
# echo -n "POSCAR."$i >> outputfiles.csv
# echo ",hey" >>outputfiles.csv
# ./coulomb.x > files2.txt
# echo "," >> files2.txt
# echo $i >> outputfiles.csv 
# echo $i",Iteration,Self,,Real,,Reciprocal" >> outputfiles.csv
# echo "Iteration,Self,,Real,,Reciprocal" >> outputfiles.csv


## code for multiple cells and same threads

# make clean > terminal.txt
# make >> terminal.txt

# cd run/
# echo "File,Reciprocal,," > parallel1.csv
# # echo "File,Self,,Reciprocal,,Real" > parallel1.csv

# for i in {1..14} 
# do
#     sed -i "s/Posfile = big\/POSCAR\.$((i-1))/Posfile = big\/POSCAR\.$i/g" input.in
#     echo -n "POSCAR"$i >> parallel1.csv
#     for j in {1..10}
#     do
#         ./coulomb.x >> parallel1.csv
#         echo " " >> parallel1.csv
#     done
#     echo " " >> parallel1.csv
#     echo " " >> parallel1.csv
# done



## code for multiple threads experiment but same cell

# cd run/
# echo "# of Threads,Self,,Reciprocal,,Real" > original_poscar15_fifty.csv
# cd ..
# for i in {2..24}
# do
#     sed -i "s/\#define NUM_THREADS $((i-1))/\#define NUM_THREADS $i/g" inc/const.h
#     make clean > output.txt
#     make >> output.txt

#     cd run/
#     echo -n $i >> original_poscar15_fifty.csv
#     for j in {1..5}
#     do
#         ./coulomb.x >> original_poscar15_fifty.csv
#         echo " " >> original_poscar15_fifty.csv
#     done
#     echo " " >> original_poscar15_fifty.csv
#     echo " " >> original_poscar15_fifty.csv
#     cd ..
# done


# # code for bench marking the calculations using lammps
# make clean > terminal.txt
# make >> terminal.txt

# cd run/
# # echo "File,Ecoul,Elong,Total" > output.csv

# for i in {1..5}
# do 
#     sed -i "s/Posfile = bench\/benchfour$((i-1))\.POSCAR/Posfile = bench\/benchfour$i\.POSCAR/g" input.in
#     # echo -n "POSCAR"$i >> output.csv
#     # for j in {1..10}
#     # do 
#         ./coulomb.x >> output1.csv
#         echo " " >> output1.csv
#     # done
#     # echo " " >> output1.csv
#     # echo " " >> output1.csv
# done

## code for multiple threads for different settings with same size and same number of atoms

for cell in {1..10}
do
    sed -i "s/Posfile = same_same\/POSCAR\.$((cell-1))/Posfile = same_same\/POSCAR\.$cell/g" run/input.in
    cd run/
    echo "# of Threads,,Reciprocal,,Real" > same_poscar$cell\file.csv
    cd ..
    for i in {2..25}
    do
        sed -i "s/\#define NUM_THREADS $((i-1))/\#define NUM_THREADS $i/g" inc/const.h
        make clean > output.txt
        make >> output.txt

        cd run/
        echo -n $i >> same_poscar$cell\file.csv
        for j in {1..5}
        do
            ./coulomb.x >> same_poscar$cell\file.csv
            echo " " >> same_poscar$cell\file.csv
        done
        echo " " >> same_poscar$cell\file.csv
        echo " " >> same_poscar$cell\file.csv
        cd ..
    done
done

# fool=1
# echo "hey" >> hey$fool\to.txt