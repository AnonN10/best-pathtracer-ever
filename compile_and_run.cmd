mkdir contrib

call :git_add_submodule https://github.com/libsdl-org/SDL, SDL, SDL2
call :git_add_submodule https://github.com/madmann91/bvh, bvh, master
call :git_add_submodule https://github.com/AnonN10/glm, glm, master

mkdir build
cd build
:: cmake .. -DSDL2_PATH="G:/C++/SDL2-2.26.2"
cmake ..
cmake --build . --config Release || pause && exit

"./Release/HelloWorld.exe"

pause
exit

:git_add_submodule
cd contrib
git submodule add %~1
cd %~2
git checkout %~3
cd ../../
EXIT /B 0