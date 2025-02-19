if [[ "$TESTOS" == "OSX" ]]; then
  cmake .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Mac-python2.cmake \
    -DCMAKE_BUILD_TYPE=Release ;
else
  if  [ ! -z ${NOIALLGATHER+x} ]; then
    cmake .. -DNOIALLGATHER=YES -DCMAKE_BUILD_TYPE=Debug ;
  elif [ ! -z ${DEBUG+x} ]; then
    cmake .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DCMAKE_BUILD_TYPE=Debug ;
  else
    cmake .. -DCMAKE_TOOLCHAIN_FILE=../Targets/Linux.cmake \
          -DCMAKE_BUILD_TYPE=Release ;
  fi
fi
