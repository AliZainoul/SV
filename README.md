# Finite Volume Saint-Venant Schema in order to APPROXIMATE the exact solution

# How to Compile:


g++ -std=c++11 -c error.cpp Vector.cpp abstract_mat_c.cpp full_mat_c.cpp mesh.cpp SV.cpp main.cpp

g++ -o m error.o Vector.o abstract_mat_c.o full_mat_c.o mesh.o SV.o main.o

./m
