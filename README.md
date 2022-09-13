# ICityTool

<img src="https://lecture11-1301936037.cos.ap-guangzhou.myqcloud.com/202209131306656.png" alt="image-2															0220913125604437" style="zoom:22%;" />	

[![PyPI version](https://img.shields.io/badge/pypi-v0.1-yellowgreen?logo=pypi&logoColor=yellow)](https://badge.fury.io/py/ICityTool) [![Python 3.6](https://img.shields.io/badge/python-3.6%7C3.7%7C3.8%7C3.9-yellowgreen?style=flat&logo=python&logoColor=yellow&color=blue)](https://badge.fury.io/py/ICityTool)[![Python 3.6](https://img.shields.io/badge/GitHub-repository-yellowgreen?style=flat&logo=github&logoColor=white&color=blue)](https://github.com/sysu-software-2022/ICityTool)

An integrate python package version of ICityRunner

## 🌟Download

```python
pip install ICityTool
```



## 👾Quick Example

```python
import ICityTool
import os
ICityTool.ICityRunner(DefenseSystem_Name="DEMO_A",
                    DefenseSystem_FilePath="./",
                    PTYFile=os.path.join("./", "DemoInput/Database/CDS.pty"),
                    PathToDatabase=os.path.join("./", "DemoInput/Database/ProteinDB"),
                    SeedPath=os.path.join("./", "DemoInput/Archaea_Cas.xlsx"),
                    NeighborhoodVicinitySize=10000,
                    PermissiveClusteringThreshold=0.3,
                    SortingOverlapThreshold=0.4,
                    SortingCoverageThresold=0.25,
                    ThreadNum="48")


```

##### I. Parameters:

1. DefenseSystem_Name: ABI, RM, TA, DND, Cas.
2. DefenseSystem_FilePath: Your working directory.
3. SeedPath: your seed **xlsx** file path
4. ThreadNum: thread number should be contingent on your **CPU core number**.





## 🧩Documentation





## 💡Reference
