git submodule update --init --recursive

mkdir build
cd build
:: cmake .. -DSDL2_PATH="G:/C++/SDL2-2.26.2"
cmake ..
cmake --build . --config Release || pause && exit

"./Release/HelloWorld.exe"

pause
exit