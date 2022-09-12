
import all
import os

if __name__ == "__main__":
    all.ICityRunner(DefenseSystem_Name="Hachiman_A",
                    DefenseSystem_FilePath="./",
                    PTYFile=os.path.join("./", "DemoInput/Database/CDS.pty"),
                    PathToDatabase=os.path.join("./", "DemoInput/Database/ProteinDB"),
                    SeedPath=os.path.join("./", "DemoInput/Archaea_Cas.xlsx"),
                    NeighborhoodVicinitySize=10000,
                    PermissiveClusteringThreshold=0.3,
                    SortingOverlapThreshold=0.4,
                    SortingCoverageThresold=0.25,
                    ThreadNum="48")
