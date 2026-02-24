The program 'MCGP' contains two input parameters: InsName and Seed.

1) InsName: the name of each instance tested.
2) Seed: the random seed.



Keep the program and the benchmark instances in the appropriate folders.

Then the job can be submitted as follows:
/************************************************\
./MCGP -i <instance_path> -s <seed>
/************************************************\

Where:
-i : path of the instance file
-s : random seed

For example:
./MCGP -i ../instances/tsplib/u159_15_5.txt -s 0



If you have any questions with this program, please contact us.

Sincere thanks and best regards!