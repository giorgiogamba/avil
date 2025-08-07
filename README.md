# avil
Audio capture and visualization in C++ for MacOS

The aim of this project is to create a program that captures the
audio coming from the pc built-in microphone and visualize
its spectrum on the screen for MacOS.

We will use standard C++ 20.

## How to run
In order to rubn this program you simply have to run `make build/avil` on your terminal, having the g++ compiler installed.
Once compiled, you can run the application by typing `./build/avil`.

## Clean
In order to clean builds, you can call `make clean` and the content of the build folder will be deleted

## Install dependencies
In order to install project dependencies, you can run `make install_deendencies` and a lib folder with the compiled libraries will be created
