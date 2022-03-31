### Python setup:
Create a local venv:
```
python -m venv venv
```
And activate it:

```
source venv/bin/activate
```

Install necessary dependencies:
```
pip install -r requirements.txt
```
### C++ lib:
Create build folder inside project:
```
mkdir build && cd build
```
#### GCC/Clang
```
cmake -DCMAKE_CXX_COMPILER=gcc .. && make
```
#### Intel
Activate icc environment variables:
```
. /opt/intel/oneapi/setvars.sh
```
Use cmake to compile lib and install it:
```
cmake -DCMAKE_CXX_COMPILER=icpc .. && make
```