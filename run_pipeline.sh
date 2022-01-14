#!/bin/bash

FILE=simulation/error.txt

if [ -f "$FILE" ]; then
  printf "$FILE exists from a previous run.\nDo you want to delete it and continue? y/n\n"
  while [ true ];
  do 
    read -p "Input Selection: " INPUT
    if [ "$INPUT" = "y" ]; then
      printf "Removing $FILE ...\n"
      rm -f "$FILE"
      printf "Start running the simulation ...\n"
      cd simulation && 
      make && 
      cd .. &&
      if [ ! -f "$FILE" ]; then
        git add RData/* simulate_all.Rout &&
        git commit -m "Push simulation results" &&
        git push
      fi
      break
    elif [ "$INPUT" = "n" ]; then
      printf "Aborting pipeline.\n"
      break
    else
      printf "Invalid selection. Input must be either 'y' or 'n'.\n"
    fi
  done
else 
  printf "Start running the simulation ...\n"
  cd simulation && 
  make && 
  cd .. &&
  if [ ! -f "$FILE" ]; then
    git add RData/* simulate_all.Rout &&
    git commit -m "Push simulation results" &&
    git push
  fi
fi
