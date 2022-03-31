### Python setup:
Create a local venv:
```
python -m venv venv
```
Install necessary dependencies:
```
pip install -r requirements.txt
```
### C++ lib compilation:
Activate icc environment variables:
```
. /opt/intel/oneapi/setvars.sh
```
Use cmake to compile lib and install it:
```
cmake -DCMAKE_CXX_COMPILER=icpc .. && make
```