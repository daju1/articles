#curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
#python3 get-pip.py
#python3 -m pip install numpy-stl

PROJECT_ROOT=${PWD}
SAGE_VER=9.1
export SAGE_ROOT=${PROJECT_ROOT}/sage-${SAGE_VER} && ${PROJECT_ROOT}/sage-${SAGE_VER}/local/bin/sage --pip install numpy-stl vtkplotlib
