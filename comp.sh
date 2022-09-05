#!/bin/bash
# Script to compile and execute a c program
gcc -ansi -Wall -Wextra -Werror -pedantic-errors kmeans.c spkmeans.c -lm -o spkmeans
