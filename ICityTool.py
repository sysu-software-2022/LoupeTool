import time
import numpy as np
import pandas as pd
import shutil
import sys
import os
import re
import subprocess
import openpyxl

from multiprocessing import Pool,set_start_method,get_start_method


def subprocess_call(Stage, CallText):

    print(Stage)
    print(CallText)
    try:
        devnull = open(os.devnull, 'wb')
        subprocess.check_call(CallText, shell=True, stdout=devnull, stderr=devnull)
        '''
        # Fetch print std
        process = subprocess.Popen(CallText, shell=True, stdout=subprocess.PIPE,stderr=subprocess.PIPE )
        while True:
            if process.poll() is None:
                line = process.stdout.readline().strip()
                print(line)
       '''
    except subprocess.CalledProcessError as e:
        print(str(e))
        print(Stage + " failed")
        sys.exit(1)


def CheckDependencies():
    DependenciesMissing = False
    if shutil.which("mmseqs") is None:
        print("mmseqs not available")
        DependenciesMissing = True

    if shutil.which("psiblast") is None:
        print("psiblast not available")
        DependenciesMissing = True

    if shutil.which("blastdbcmd") is None:
        print("blastdbcmd not available")
        DependenciesMissing = True

    if shutil.which("python") is None:
        print("python not available")
        DependenciesMissing = True

    if shutil.which("muscle") is None:
        print("muscle not available")
        DependenciesMissing = True

    if DependenciesMissing:
        print("required software missing")
        sys.exit(1)

    return


# step 5
def HashSeedExtract(cds_path="Database/archaea_cds.pty", seed_path="Archaea_RM.xlsx", seeds_name=os.path.join("./" + "RM_A"+"_OUTPUT", "Seeds_" + "RM_A" + ".tsv")):
    # Hyper-parameters
    global i
    i = 0
    global cnt
    cnt = 0
    global read_time_e
    read_time_e = 0
    global read_time_s
    read_time_s = 0
    global t1
    t1 = 0
    global t2
    t2 = 0

    # O(n) Create Hash table
    def create_hash_table(cds_find):
        hash_table = dict()
        for idx, row in enumerate(cds_find.itertuples()):
            if getattr(row, 'cds_accession') not in hash_table.keys():
                hash_table.update({getattr(row, 'cds_accession'): idx + 1})
        return hash_table


    def SeedExtract(seed_path, cds_path):
        global i, read_time_s, read_time_e, t1, t2
        info_alloc = pd.DataFrame()
        read_time_s = time.time()
        cds_file = pd.read_table(cds_path, names=['LociID', 'ORFStart_ORFStop',
                                                'Strand', 'OrganismID', 'ContigID', 'product_accession'])

        cds_file[['Start', 'Stop']] = cds_file["ORFStart_ORFStop"].str.split(':', n=1, expand=True)
        cds_file['Stop'] = cds_file['Stop'].str.replace('>', '')
        cds_file[['OrganismID', 'assembly_accession']] = cds_file['OrganismID'].str.rsplit('-', n=1, expand=True)
        cds_file = cds_file.drop(labels='ORFStart_ORFStop', axis=1)
        cds_file = cds_file.drop(labels='OrganismID', axis=1)

        cds_accession = cds_file.assembly_accession
        cds_end = cds_file.Stop
        cds_find = pd.concat([cds_accession, cds_end], axis=1); cds_find.columns = ['cds_accession', 'cds_end']

        hash_table = create_hash_table(cds_find)
        hash_table_values = list(hash_table.values())
        print("----------------------------HASH TABLE----------------------------",flush=True)
        print(hash_table, flush=True)
        print("HASH_TABLE LENGTH = ", len(hash_table), flush=True)

        cds_find = cds_find.to_numpy()

        seed_file = pd.read_excel(seed_path)
        seed_assembly_accession = seed_file.assembly_accession.str.replace("GCA", "GCF")
        seed_end = seed_file.end
        seed_find = pd.concat([seed_assembly_accession, seed_end], axis=1)

        read_time_e = time.time()

        print("------------------------Processing Starting Now------------------------", flush=True)
        t1 = time.time()
        for id, end in zip(seed_find.assembly_accession, seed_find.end):
            if id not in hash_table.keys():
                continue
            if hash_table_values.index(hash_table[id]) + 1 < len(hash_table_values):
                l2 = hash_table_values[hash_table_values.index(hash_table[id]) + 1]
            else:
                l2 = len(cds_file)
            l1 = int(hash_table[id])  # Starting Point
            length = l2 - l1
            print(id, ' ---------> Search length =', length, '[Start =', l1, ', Stop =', l2, ']', flush=True)
            accession_list = cds_find[int(hash_table[id]): int(hash_table[id]) + length + 1]

            for idx, item in enumerate(accession_list):
                if (item == np.array([id, str(end)])).all():
                    print('No.', i, ' Seed Found: ', id, str(end), flush=True)
                    # Process the corresponding items
                    data = cds_file.loc[idx + l1, ['assembly_accession',
                                            'LociID', 'product_accession', 'ContigID', 'Start', 'Stop']]
                    data = data.to_dict()
                    data = pd.DataFrame(data, index=[i]); i += 1
                    info_alloc = pd.concat([info_alloc, data])
            t2 = time.time()
        info_alloc.drop_duplicates(inplace=True)
        print("------------------------Processing Finished------------------------",flush=True)
        print(info_alloc,flush=True)
        print('----------------------------', i, 'Seeds in total----------------------------',flush=True)
        info_alloc.to_csv(seeds_name, index=False, sep='\t', header=False, )


    # Hit it!
    SeedExtract(seed_path, cds_path)
    print("---------------------------------Performance Attributes---------------------------------")
    print('File saved as ' + seeds_name)
    print('Reading files and creating hash table used time = %.3fs' %(read_time_e - read_time_s))
    print('Searching used time = %.3fs' %(t2 - t1))
    print("Total used time = %.3fs" %(read_time_e - read_time_s + t2 - t1))
    print("Seed Extraction Done!!")


# step 6
def FindNeighborhood(PTYDataFileName="Database/archaea_cds.pty", SeedFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "Seeds_" + "RM_A" + ".tsv"), ResultFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "Vicinity_"+"RM_A"+".tsv"), SeedIslandOffset=10000):
    SeedIslandOffset = int(SeedIslandOffset)
    def LoadSeedDirt(SeedFileName, NoID =False):
        Seeds = dict()

        ContigID = 3
        Start = 4
        End = 5
        Offset = 0

        return AddSeedDictFields(SeedFileName, Seeds, ContigID, Start, End, Offset, NoID)

    def AddSeedDictFields(SeedFileName, Seeds, ContigID, Start, End, Offset, NoID=False):
        with open(SeedFileName, 'r') as File:
            for Line in File:
                LineValues = Line.split('\t')
                LineValues[-1] = LineValues[-1][:-1]

                if not NoID:
                    ID = LineValues[2]
                else:
                    ID = 'Seed'

                ORF_start, ORF_end = minmax(int(LineValues[Start]), int(LineValues[End]))
                if LineValues[ContigID] in Seeds:
                    Seeds[LineValues[ContigID]].append([max(ORF_start - Offset, 0), ORF_end + Offset, ID, LineValues])
                else:
                    Seeds[LineValues[ContigID]] = [[max(ORF_start - Offset, 0), ORF_end + Offset, ID, LineValues]]
                
            for Key in Seeds:
                Seeds[Key] = sorted(Seeds[Key], key= lambda e:e[0])
        return Seeds

    def minmax(X, Y):
        return min(X, Y), max(X, Y)

    # CDS PTY tsv file columns
    global PTY_LocusTag
    global PTY_Coordinates
    global PTY_Strand
    global PTY_SpecieName
    global PTY_ContigId
    global PTY_AccessionNo

    PTY_LocusTag = 0
    PTY_Coordinates = 1
    PTY_Strand = 2
    PTY_SpecieName = 3
    PTY_ContigId = 4
    PTY_AccessionNo = 5

    def GetSeed(ContigID, Start, End, Seeds):
        if not ContigID in Seeds:
            return []
        
        for ORF in Seeds[ContigID]:
            if(Start <= ORF[1]) and (End >= ORF[0]):
                return ORF

        return []

    def IsInSeeds(ContigID, Start, End, Seeds):
        InSeed = GetSeed(ContigID, Start, End, Seeds)
        if len(InSeed) == 0:
            return False
        return True

    def WriteSeedIslands(SeedIslandPTYLines, ResultFile):

        if len(SeedIslandPTYLines) > 0:
            ResultFile.write('===\n')
            for LineValues in SeedIslandPTYLines:
                ResultFile.write(
                LineValues[0] + '\t' + 
                LineValues[1].replace(':','..') + '\t' +
                LineValues[2] + '\t' + 
                LineValues[3] + '\t' +
                LineValues[4] + '\n')


    def WriteIslands(ContigLines, ResultFile, ContigID, Seeds, SeedIslandOffset):
        if ContigID not in Seeds:
            return
        if len(Seeds[ContigID]) == 0:
            return 
        
        
        ContigLines = sorted(ContigLines, key= lambda e:e[0])

        SeedIslandPTYLines = []
        for Line in ContigLines:
            LineValues = Line[1]

            Coordinates = LineValues[PTY_Coordinates].split(":")
            
            
            try:
                Start = int(Coordinates[0])
                End = int(Coordinates[1])
            except ValueError:
                continue

            if IsInSeeds(ContigID, Start - SeedIslandOffset, End + SeedIslandOffset, Seeds):
                SeedIslandPTYLines.append(
                    [LineValues[PTY_AccessionNo], LineValues[PTY_Coordinates], LineValues[PTY_Strand],
                    LineValues[PTY_SpecieName], LineValues[PTY_ContigId]])
            else:
                if len(SeedIslandPTYLines) > 0:
                    WriteSeedIslands(SeedIslandPTYLines, ResultFile)

                SeedIslandPTYLines = []
        
        if len(SeedIslandPTYLines) > 0:
            WriteSeedIslands(SeedIslandPTYLines, ResultFile)

    Seeds = LoadSeedDirt(SeedFileName)

    ContigID = ''
    ContigLines = []
    with open(ResultFileName, 'w') as ResultFile:
        with open(PTYDataFileName,'r') as PTYData:
            for Line in PTYData:
                Line = Line[:-1]
                LineValues = Line.split('\t')

                if ContigID != LineValues[PTY_ContigId]:
                    WriteIslands(ContigLines, ResultFile, ContigID, Seeds, SeedIslandOffset)
                    ContigLines =[]

                ContigID = LineValues[PTY_ContigId]
                Coordinates = LineValues[PTY_Coordinates].split(':')
                Coordinates[0] = re.findall(r'\d+', Coordinates[0])[0]
                Start = int(Coordinates[0])


                ContigLines.append([Start, LineValues])

            WriteIslands(ContigLines, ResultFile, ContigID, Seeds, SeedIslandOffset)


# step 7
def CollectingProteinIDs(VicinityFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "Vicinity_"+"RM_A"+".tsv"),
                        VicinityIDsFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "VicinityIDs_"+"RM_A"+".lst")):
    subprocess_call("Step 7: Collecting protein IDs", "grep -v \"===\" " + VicinityFileName + " | cut -f1 | sort -u > " +
                VicinityIDsFileName)


# step 8
def FetchingProteinSequences(PathToDatabase="Database/ArchaeaProt",
                            VicinityIDsFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "VicinityIDs_"+"RM_A"+".lst"),
                            VicinityFASTAFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "Vicinity_"+"RM_A"+".faa")):
    subprocess_call("Step 8: Fetching protein sequences", "blastdbcmd -db " + PathToDatabase +
                " -entry_batch " + VicinityIDsFileName +
                " -long_seqids > " + VicinityFASTAFileName)


# step 9
def ClusteringProteinSeqiences(VicinityFASTAFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "Vicinity_"+"RM_A"+".faa"),
                            PermissiveClusteringThreshold=0.3,
                            VicinityClustersFileName="VicinityPermissiveClustsLinear" + "RM_A" + ".tsv"):

    subprocess_call("Step 9: Clustering protein seqiences", "bash RunClust.sh " + VicinityFASTAFileName + " " +
                str(PermissiveClusteringThreshold) + " " +
                VicinityClustersFileName)


# step 10
def generator(ClusterNo, Sequence, Database, ClustersFolderName):
    ClusterIDs = Sequence[:-1].split("\t")[1].split(" ")
    writer(ClusterIDs, ClusterNo, Database, ClustersFolderName)
    return Sequence

def writer(ClusterIDs, ClusterNo, Database, ClustersFolderName):
    global cpu_cnt
    ClusterProfileFileName = "CLUSTER_" + str(ClusterNo + 1) + ".ali"
    
    """
    Critical Section
    """

    TmpIDsFileName = "Tmp_IDs" + "_" + str(ClusterNo) + ".lst"
    TmpFASTAFileName = "Tmp_FASTA" + "_" + str(ClusterNo) + ".faa"

    with open(TmpIDsFileName, "w") as IDsFile:
        IDsFile.write("\n".join(ClusterIDs))

    subprocess.call("blastdbcmd -db " + Database + " -entry_batch " + TmpIDsFileName + " -long_seqids > " + TmpFASTAFileName
                    , shell=True)

    subprocess.call("muscle -align " + TmpFASTAFileName + " -output " + ClustersFolderName + "/" + ClusterProfileFileName
                    , shell=True)

    
    subprocess.call("rm " + TmpIDsFileName, shell=True)

    subprocess.call("rm " + TmpFASTAFileName, shell=True)


def MakeProfiles(ClustersFileName="VicinityPermissiveClustsLinear" + "RM_A" + ".tsv", ClustersFolderName=os.path.join("./" + "RM_A"+"_OUTPUT", "CLUSTERS_"+"RM_A" + "/"), Database="RM_A"):

    global cpu_cnt
    cpu_cnt = os.cpu_count()
    ClusterNo = 0

    if not os.path.exists(ClustersFolderName):
        subprocess.call("mkdir " + ClustersFolderName, shell = True)


    def multi_process():
        global cpu_cnt
        if not os.path.getsize(ClustersFileName):
            print("Clusters File empty")
            exit()

        with open(ClustersFileName, 'r') as File, Pool(cpu_cnt) as pool:
            Sequences = File.readlines(); Length = len(Sequences)
            for ClusterNo, Sequence in enumerate(Sequences):
                pool.apply_async(generator, (ClusterNo, Sequence, Database, ClustersFolderName))
            pool.close()
            pool.join()


    if __name__ == "all":
        t1 = time.time()
        cnt = 0

        multi_process()

        t2 = time.time()

        print("Work Done !")
        print("Used time: ", t2 - t1, "s")


# step 11
def RunPSIBLAST(ClustersFolderName=os.path.join("./" + "RM_A"+"_OUTPUT", "CLUSTERS_"+"RM_A" + "/"), Database="Database/ArchaeaProt", ThreadNum="48"):

    def GetNumberOfSequences(FASTAFileName):
        Count = 0
        for Line in open(FASTAFileName):
            if ">" in Line:
                Count += 1

        return Count



    global Evalue
    Evalue = 0.000001

    for FileName in os.listdir(ClustersFolderName):
        if not FileName.endswith("ali"):
            continue

        BLASTHitsFileName = ClustersFolderName + "/" + os.path.splitext(FileName)[0] + ".hits"

        if GetNumberOfSequences(ClustersFolderName + "/" + FileName) == 1:
            Query = "-query"
        else:
            Query = "-in_msa"
        Query += " " + ClustersFolderName + "/" + FileName

        # to parallelize
        print("Processing " + FileName)
        subprocess.call(
            "psiblast -db " + Database + " " +
            "-outfmt \"7 qseqid sseqid slen sstart send evalue qseq sseq qstart qend score\"" +
            " -seg no -evalue " + str(Evalue) +
            " -dbsize 20000000 -max_target_seqs 10000 -comp_based_stats no " +
            "-num_threads " + ThreadNum + " " +
            Query + " > " + BLASTHitsFileName, shell = True)


# step 12
def SortBLASTHitsInMemory(ClustersHitsFolder=os.path.join("./" + "RM_A"+"_OUTPUT", "CLUSTERS_"+"RM_A" + "/"),
                        SortedHitsFolder=os.path.join("./" + "RM_A"+"_OUTPUT", "CLUSTERS_"+"RM_A", "Sorted/"),
                        PTYFileName="Database/archaea_cds.pty",
                        VicinityProteinIDsFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "VicinityIDs_"+"RM_A"+".lst"),
                        SeedsFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "Seeds_" + "RM_A" + ".tsv"),
                        VicinityFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "Vicinity_"+"RM_A"+".tsv"),
                        MIN_OVERLAP_FILTER=0.4,
                        MIN_COVERAGE_FILTER=0.25):
    MIN_OVERLAP_FILTER = float(MIN_OVERLAP_FILTER)
    MIN_COVERAGE_FILTER = float(MIN_COVERAGE_FILTER)
    # Reading arguments from command line
    global HIT_Score
    global HIT_Cluster
    global HIT_ALIGNMENT_START
    global HIT_ALIGNMENT_STOP
    global HIT_ALIGNMENT_Sequence
    global HIT_Protein_ID

    global INFO_CONTIG
    global INFO_BAITED
    global INFO_ORF_START
    global INFO_ORF_STOP
    global INFO_ORF_DISTANCE_TO_BAIT

    HIT_Score = 1
    HIT_Cluster = 5
    HIT_ALIGNMENT_START = 2
    HIT_ALIGNMENT_STOP = 3
    HIT_ALIGNMENT_Sequence = 4
    HIT_Protein_ID = 0

    INFO_CONTIG = 0
    INFO_BAITED = 1
    INFO_ORF_START = 2
    INFO_ORF_STOP = 3
    INFO_ORF_DISTANCE_TO_BAIT = 4


    def AddHitToDict(TargetHitsDict, Hit):
        BetterHitExists = False

        for DictHit in TargetHitsDict[Hit[HIT_Protein_ID]]:
            InterCoverage = min(Hit[HIT_ALIGNMENT_STOP],
                                DictHit[HIT_ALIGNMENT_STOP] - max(Hit[HIT_ALIGNMENT_START], DictHit[HIT_ALIGNMENT_START]))
            MinLength = min(Hit[HIT_ALIGNMENT_STOP] - Hit[HIT_ALIGNMENT_START],
                            DictHit[HIT_ALIGNMENT_STOP] - DictHit[HIT_ALIGNMENT_START])

            if InterCoverage / MinLength > MIN_COVERAGE_FILTER:
                if Hit[HIT_Score] > DictHit[HIT_Score]:
                    TargetHitsDict[Hit[HIT_Protein_ID]].remove(DictHit)
                else:
                    BetterHitExists = True
                    break

        if not BetterHitExists:
            TargetHitsDict[Hit[HIT_Protein_ID]].append(Hit)

        return TargetHitsDict



    def AddTargetHitsValues(TargetHitsFileName, TargetHitsDict, ProteinInfoDict):
        BLAST_FORMAT_ALIGNMENT_TARGET_PROTEIN_ID = 1
        BLAST_FORMAT_ALIGNMENT_START = 8
        BLAST_FORMAT_ALIGNMENT_STOP = 9
        BLAST_FORMAT_ALIGNMENT_SCORE = 10
        BLAST_FORMAT_ALIGNMENT_SEQUENCE = 7

        TargetValues = []
        BestScoreLine = []
        MaxScore = 0
        ClusterId = os.path.splitext(os.path.basename(TargetHitsFileName))[0]

        with open(TargetHitsFileName, "r") as TargetHitsFile:
            for Line in TargetHitsFile:
                if Line[0] == "#":
                    continue

                LineValues = Line.split("\t")
                LineValues[-1] = LineValues[-1][:-1]

                # Alignment Data
                Start = int(LineValues[BLAST_FORMAT_ALIGNMENT_START])
                Stop = int(LineValues[BLAST_FORMAT_ALIGNMENT_STOP])
                Score = int(LineValues[BLAST_FORMAT_ALIGNMENT_SCORE])
                if "|" in LineValues[BLAST_FORMAT_ALIGNMENT_TARGET_PROTEIN_ID]:
                    ProteinID = LineValues[BLAST_FORMAT_ALIGNMENT_TARGET_PROTEIN_ID].split("|")[
                        1]  # Accession|ID| format expected
                else:
                    ProteinID = LineValues[BLAST_FORMAT_ALIGNMENT_TARGET_PROTEIN_ID]

                HitLine = [ProteinID, Score, Start, Stop, LineValues[BLAST_FORMAT_ALIGNMENT_SEQUENCE], ClusterId]
                HitLine.extend(ProteinInfoDict[ProteinID])

                if Score > MaxScore:
                    BestScoreLine = HitLine
                    MaxScore = Score

                TargetValues.append(HitLine)

        if len(TargetValues) > 0:  # no hits
            # max length of alignment is taken as selflength to filter by coverage
            # this value is used as representative length of alignment
            BestScoreLength = len(BestScoreLine[HIT_ALIGNMENT_Sequence].replace("-", ""))

            for Hit in TargetValues:
                if len(Hit[HIT_ALIGNMENT_Sequence].replace("-",
                                                        "")) / BestScoreLength >= MIN_OVERLAP_FILTER:  # filtering by coverage
                    if Hit[HIT_Protein_ID] in TargetHitsDict:
                        TargetHitsDict = AddHitToDict(TargetHitsDict, Hit)
                    else:
                        TargetHitsDict[Hit[HIT_Protein_ID]] = [Hit]

        return ClusterId

    def GetDistToSeedDict(AnnotationFileName, SeedsDict):
        LOCI_PROTEIN_ID = 0
        LOCI_CONTIG_ID = 4
        LOCI_COORDINATES = 1

        DistToSeed = dict()

        Accessions = dict()
        IslandSeeds = []
        for Line in open(AnnotationFileName, "r"):
            LineValues = Line[:-1].split("\t")

            if LineValues[LOCI_PROTEIN_ID] == "===":
                if len(Accessions) > 0:
                    if len(IslandSeeds) == 0:
                        raise "NoSeedFound"
                    else:
                        for Accession in Accessions:
                            DistToSeed[Accession] = min([abs(x - Accessions[Accession]) for x in IslandSeeds])

                Accessions = dict()
                IslandSeeds = []
            else:
                Accessions[LineValues[LOCI_PROTEIN_ID]] = len(Accessions) + 1

                StartStop = LineValues[LOCI_COORDINATES].split("..")
                if IsInSeeds(LineValues[LOCI_CONTIG_ID], int(StartStop[0]), int(StartStop[1]), SeedsDict):
                    IslandSeeds.append(len(Accessions))

        if len(Accessions) > 0:
            if len(IslandSeeds) == 0:
                raise "NoSeedFound"
            else:
                for Accession in Accessions:
                    DistToSeed[Accession] = min([abs(x - Accessions[Accession]) for x in IslandSeeds])

        return DistToSeed

    def LoadProteinInfoDict(PTYFileName, VicinityProteinIDs, LociDists):
        PTY_PROTEIN_ID = 5
        PTY_CONTIG_ID = 4
        PTY_COORDINATES = 1

        ProteinInfoDict = dict()

        for Line in open(PTYFileName):
            LineValues = Line[:-1].split("\t")

            StartStop = LineValues[PTY_COORDINATES].split(":")

            if LineValues[PTY_PROTEIN_ID] in VicinityProteinIDs:
                InVicinity = "1"
            else:
                InVicinity = "0"

            if LineValues[PTY_PROTEIN_ID] in LociDists:
                Dist = LociDists[LineValues[PTY_PROTEIN_ID]]
            else:
                Dist = 10000

            # in assumption that Protein IDs are unique
            ProteinInfoDict[LineValues[PTY_PROTEIN_ID]] = [LineValues[PTY_CONTIG_ID], InVicinity,
                                                        StartStop[0], StartStop[1], Dist]

        return ProteinInfoDict

    def LoadList(FileName):
        ListValues = []
        for Line in open(FileName):
            ListValues.append(Line[:-1])
        return ListValues

    def LoadSeedsDict(SeedsFileName):
        Seeds = dict()
        ContigField = 3
        StartField = 4
        EndField = 5
        Offset = 0

        return AddSeedsDictByFields(SeedsFileName, Seeds, ContigField, StartField, EndField, Offset)

    def AddSeedsDictByFields(FileName, SeedDict, ContigField, StartField, EndField, Offset):
        with open(FileName, "r") as File:
            for Line in File:
                LineValues = Line.split("\t")
                LineValues[-1] = LineValues[-1][:-1]
                ID = LineValues[0]

                Start, End = minmax(int(LineValues[StartField]), int(LineValues[EndField]))
                if LineValues[ContigField] in SeedDict:
                    SeedDict[LineValues[ContigField]].append([max(Start - Offset, 0), End + Offset, ID, LineValues])
                else:
                    SeedDict[LineValues[ContigField]] = [[max(Start - Offset, 0), End + Offset, ID, LineValues]]

        for Key in SeedDict:
            SeedDict[Key] = sorted(SeedDict[Key], key = lambda e: e[0])

        return SeedDict

    def minmax(X, Y):
        return min(X, Y), max(X, Y)

    def IsInSeeds(ContiAccessionD, Start, Stop, SeedsDict):
        Seed = GetSeed(ContiAccessionD, Start, Stop, SeedsDict)
        if len(Seed) == 0:
            return False
        return True

    def GetSeed(ContiAccessionD, Start, Stop, SeedsDict):
        if not ContiAccessionD in SeedsDict:
            return []

        for ORF in SeedsDict[ContiAccessionD]:
            if (Start <= ORF[1]) and (Stop >= ORF[0]):
                return ORF

        return []

    VicinityProteinIDs = set(LoadList(VicinityProteinIDsFileName))
    Seeds = LoadSeedsDict(SeedsFileName)
    LociDists = GetDistToSeedDict(VicinityFileName, Seeds)
    ProteinInfoDict = LoadProteinInfoDict(PTYFileName, VicinityProteinIDs, LociDists)

    if not os.path.exists(SortedHitsFolder):
        os.makedirs(SortedHitsFolder)

    ProteinHitsDict = dict()
    Clusters = []
    for FileName in os.listdir(ClustersHitsFolder):
        if not FileName.endswith(".hits"):
            continue

        ClusterId = AddTargetHitsValues(ClustersHitsFolder + FileName, ProteinHitsDict, ProteinInfoDict)
        Clusters.append(ClusterId)

    print("Loading BLAST hits complete.")

    ClusterHitsDict = dict()
    for ProteinID in ProteinHitsDict:
        for Hit in ProteinHitsDict[ProteinID]:
            if Hit[HIT_Cluster] in ClusterHitsDict:
                ClusterHitsDict[Hit[HIT_Cluster]].append(Hit)
            else:
                ClusterHitsDict[Hit[HIT_Cluster]] = [Hit]

    print("Moving sorting hits by clusters complete.")

    for ClusterId in Clusters:
        if ClusterId in ClusterHitsDict:
            with open(SortedHitsFolder + str(ClusterId) + ".hits_sorted", "w") as SortedFile:
                for Hit in ClusterHitsDict[ClusterId]:
                    ResLine = "\t".join([str(x) for x in Hit])

                    SortedFile.write(ResLine + "\n")
        else: # empty file
            with open(SortedHitsFolder + str(ClusterId) + ".hits_sorted", "w") as SortedFile:
                SortedFile.write("")


#step 13
def CalculatingICITYMetric(SortedBLASTHitsFolder=os.path.join("./" + "RM_A"+"_OUTPUT", "CLUSTERS_"+"RM_A", "Sorted/"),
                        PathToDatabase="Database/ArchaeaProt",
                        VicinityClustersFileName="VicinityPermissiveClustsLinear" + "RM_A" + ".tsv",
                        ICITYFileName=os.path.join("./" + "RM_A"+"_OUTPUT", "Relevance_"+"RM_A"+".tsv"),
                        ThreadNum="48"):
    subprocess_call("Step 13: Calculating ICITY metric", "bash Cal.sh " +
                    SortedBLASTHitsFolder + " " +
                    PathToDatabase + " " +
                    VicinityClustersFileName + " " +
                    ICITYFileName + " " +
                    ThreadNum
                    )


# step 14
def SortRelevance(DefenseSystem_Name="RM_A", DefenseSystem_FilePath="./"):
    # Relevance.tsv
    CLUSTER_NAME = 0
    NUMBER_IN_VICINTY = 1
    NUMBER_IN_ENTIRE = 2
    DISTANCE_TO_SEED = 3
    DS3P = 4




    DS3P_threshold = 0.0


    PROTEIN_ACCESSION = 0

    SEED_ACCESSION = 2


    SeedAccession_Filename = DefenseSystem_FilePath+ DefenseSystem_Name + '_OUTPUT/Seeds_' + DefenseSystem_Name + '.tsv'
    CLUSTERS_Filepath = DefenseSystem_FilePath + DefenseSystem_Name + '_OUTPUT/CLUSTERS_' + DefenseSystem_Name
    Vicinity_Filename = (os.path.split(DefenseSystem_FilePath+ DefenseSystem_Name + '_OUTPUT/Vicinity_'+ DefenseSystem_Name + '.faa')[1])
    Relevance_Filenpath = DefenseSystem_FilePath + DefenseSystem_Name + '_OUTPUT/Relevance_' + DefenseSystem_Name + '.tsv'
    Relevance_CategoryName = 'Relevance_Sorted_'+ DefenseSystem_Name + '_Category.csv'
    CLUSTERHitsSorted_path = DefenseSystem_FilePath + DefenseSystem_Name + '_OUTPUT/CLUSTERS_' + DefenseSystem_Name + '/Sorted/'


    DefenseSystem_Function_1 = 'abortive infection'
    DefenseSystem_Function_2 = 'restriction-modification'
    DefenseSystem_Function_3 = 'toxin'
    DefenseSystem_Function_4 = 'argonaute'
    DefenseSystem_Function_5 = 'Dnd'
    Dir = "ACCESSION_A"

    if not os.path.exists(Dir):
        os.mkdir(Dir)


    subprocess.call("cat " + CLUSTERS_Filepath + "/CLUSTER_*.ali  | grep ref | awk -F \"|\" '{print $2, \"|\", $NF}'" + " > ACCESSION_A/ACCESSION_"+ DefenseSystem_Name+'.txt',
    shell= True)
    subprocess.call('cat ACCESSION_A/ACCESSION_'+ DefenseSystem_Name + '.txt | grep \''+ DefenseSystem_Function_1 +'\' | awk \'{print $0\" |  1\" }\' > ACCESSION_A/ACCESSION_ONLY_'+ DefenseSystem_Name +'.txt',
    shell= True)
    subprocess.call('cat ACCESSION_A/ACCESSION_' + DefenseSystem_Name + '.txt | grep -E \"'+ DefenseSystem_Function_2 + '|CRISPR|'+ DefenseSystem_Function_3 +'|'+DefenseSystem_Function_4+'|'+DefenseSystem_Function_5+'\" | awk \'{print $0  \" |  2\" }\' > ACCESSION_A/ACCESSION_Other_DefenseGene_'+ DefenseSystem_Name +'.txt',
    shell= True)
    subprocess.call('cat ACCESSION_A/ACCESSION_'+ DefenseSystem_Name +'.txt | grep -v \"toxin\|CRISPR\|abortive infection\|restriction-modification\|hypothetical\|argonaute\Dnd" | awk \'{print $0  \" |  3\" }\' > ACCESSION_A/ACCESSION_HouseKeepingGene_'+ DefenseSystem_Name + '.txt',
    shell= True)
    subprocess.call('cat ACCESSION_A/ACCESSION_'+ DefenseSystem_Name + '.txt | grep \'hypothetical\' | awk \'{print $0  \" |  4\" }\' > ACCESSION_A/ACCESSION_hypothetical_'+ DefenseSystem_Name + '.txt',
    shell= True)




    Accession_All_FileName = 'ACCESSION_A/ACCESSION_' + DefenseSystem_Name + '.txt'
    Accession_Category_Target_FileName = 'ACCESSION_A/ACCESSION_ONLY_'+ DefenseSystem_Name + '.txt'
    Accession_Category_Hypothetical_FileName = 'ACCESSION_A/ACCESSION_hypothetical_'+ DefenseSystem_Name + '.txt'
    ACCESSION_HouseKeepingGene_FileName = 'ACCESSION_A/ACCESSION_HouseKeepingGene_'+ DefenseSystem_Name + '.txt'
    ACCESSION_Other_DefenseGene_FileName = 'ACCESSION_A/ACCESSION_Other_DefenseGene_'+ DefenseSystem_Name + '.txt'

    with open(Relevance_Filenpath,'r') as F:
        Relevance = dict()
        for Line in F:
            LineValue = Line[:-1].split('\t')
            Relevance[LineValue[CLUSTER_NAME]] = [LineValue[CLUSTER_NAME], LineValue[NUMBER_IN_VICINTY], LineValue[NUMBER_IN_ENTIRE], LineValue[DISTANCE_TO_SEED], LineValue[DS3P]]

    Relevance_Sorted_Target = dict()
    for r in Relevance:
        if int(Relevance[r][NUMBER_IN_ENTIRE]) >= 1:
            if float(Relevance[r][DS3P]) > DS3P_threshold:
                Relevance_Sorted_Target[r] = Relevance[r]
        else:
            continue


    for n in Relevance_Sorted_Target:
        CLUSTERHitsSorted_Filename = CLUSTERHitsSorted_path + n + '.hits_sorted'
        with open(CLUSTERHitsSorted_Filename,'r') as F:
            ProteinAccession = []
            for Line in F:
                LineValue = Line[:-1].split('\t')
                ProteinAccession.append(LineValue[PROTEIN_ACCESSION])
        
        Relevance_Sorted_Target[n] = [Relevance_Sorted_Target[n], ProteinAccession]



    Accession_Category_Target = dict()
    Accession_Category_Hypothetical = dict()
    Accession_Category_OtherDefense = dict()
    Accession_Category_HouseKeeping = dict()

    with open(Accession_Category_Target_FileName,'r') as T:
        for Line in T:
            LineValue = Line[:-1].split(' |  ')
            Accession_Category_Target[LineValue[0]] = [LineValue[1], LineValue[2]]

    with open(ACCESSION_Other_DefenseGene_FileName,'r') as T:
        for Line in T:
            LineValue = Line[:-1].split(' |  ')
            Accession_Category_HouseKeeping[LineValue[0]] = [LineValue[1], LineValue[2]]

    with open(ACCESSION_HouseKeepingGene_FileName,'r') as T:
        for Line in T:
            LineValue = Line[:-1].split(' |  ')
            Accession_Category_OtherDefense[LineValue[0]] = [LineValue[1], LineValue[2]]

    with open(Accession_Category_Hypothetical_FileName,'r') as T:
        for Line in T:
            LineValue = Line[:-1].split(' |  ')
            Accession_Category_Hypothetical[LineValue[0]] = [LineValue[1], LineValue[2]]

    Relevance_Sorted_Target_Category = dict()
    for id in Relevance_Sorted_Target:
        Accession_Category= {}
        for ac in Relevance_Sorted_Target[id][1]:
            if ac in Accession_Category_Target:
                Accession_Category[ac] = 1
                continue
            elif ac in Accession_Category_OtherDefense:
                Accession_Category[ac] = 2
                continue
            elif ac in Accession_Category_HouseKeeping:
                Accession_Category[ac] = 3
                continue
            elif ac in Accession_Category_Hypothetical:
                Accession_Category[ac] = 4
                continue
            else:
                Accession_Category[ac] = 5
        Cluster_Category_Name = min(Accession_Category, key= Accession_Category.get)
        Cluster_Category_Value = min(Accession_Category.values())
        if Cluster_Category_Value == 1:
            Relevance_Sorted_Target_Category[id] = [Relevance_Sorted_Target[id][0], Accession_Category_Target[Cluster_Category_Name][1]]
        elif Cluster_Category_Value == 2:
            Relevance_Sorted_Target_Category[id] = [Relevance_Sorted_Target[id][0], Accession_Category_OtherDefense[Cluster_Category_Name][1]]
        elif Cluster_Category_Value == 3:
            Relevance_Sorted_Target_Category[id] = [Relevance_Sorted_Target[id][0], Accession_Category_HouseKeeping[Cluster_Category_Name][1]]
        elif Cluster_Category_Value == 4:
            Relevance_Sorted_Target_Category[id] = [Relevance_Sorted_Target[id][0], Accession_Category_Hypothetical[Cluster_Category_Name][1]]


    def CalConservation(Accession_Taxonomy):
            count = 0
            if not 'None' in Accession_Taxonomy:
                for key in Accession_Taxonomy.keys():
                    count = count + 1/len(Accession_Taxonomy[key])
                return count/len(Accession_Taxonomy)
            else:
                for key in Accession_Taxonomy.keys():
                    if key != 'None':
                        count = count + 1/len(Accession_Taxonomy[key])
                return count/len(Accession_Taxonomy)


    Relevance_Taxonomy = dict()
    Relevance_Sorted_Target_Taxonomy = dict()

    with open(Accession_All_FileName,'r') as T:
        for Line in T:
            LineValue = Line[:-1].split(' |  ')
            TaxonomyID = LineValue[1].split('[')[1][:-1]
            TaxonomyID = TaxonomyID.split(' ',1)
            Relevance_Taxonomy[LineValue[0]] = TaxonomyID

    for id in Relevance_Sorted_Target:
        Accession_Taxonomy_genus = {}
        Accession_Taxonomy_species = {}

        for ac in Relevance_Sorted_Target[id][1]:
            if ac in Relevance_Taxonomy:
                Accession_Taxonomy_genus[ac] = Relevance_Taxonomy[ac][0]
                try:
                    Accession_Taxonomy_species[ac] = Relevance_Taxonomy[ac][1]
                except IndexError:
                    Accession_Taxonomy_species[ac] = 'None'

            flipped_genus = {}
            flipped_species = {}

            for key,value in Accession_Taxonomy_genus.items():
                if value not in flipped_genus:
                    flipped_genus[value] = list()
                    flipped_genus[value].append(key)
                else:
                    flipped_genus[value].append(key)

            for key,value in Accession_Taxonomy_species.items():   
                if value not in flipped_species:
                    flipped_species[value] = list()
                    flipped_species[value].append(key)
                else:
                    flipped_species[value].append(key)
        if flipped_genus == {}:
            continue
        else:
            Relevance_Sorted_Target_Taxonomy[id] = [CalConservation(flipped_genus),CalConservation(flipped_species)]

    with open(Relevance_CategoryName,'w') as f:
        for i in Relevance_Sorted_Target_Category:
            for element in Relevance_Sorted_Target_Category[i][0]:
                f.write(element)
                f.write(',')
            f.write(Relevance_Sorted_Target_Category[i][1])
            for n in Relevance_Sorted_Target_Taxonomy[i]:
                f.write(',')
                f.write(str(n))
                
            f.write("\n")



def ICityRunner(DefenseSystem_Name,
                DefenseSystem_FilePath,
                PTYFile,
                PathToDatabase,
                SeedPath,
                NeighborhoodVicinitySize,
                PermissiveClusteringThreshold,
                SortingOverlapThreshold,
                SortingCoverageThresold,
                ThreadNum):
    CheckDependencies()
    Dir = DefenseSystem_Name + "_OUTPUT"
    if not os.path.exists(Dir):
        os.mkdir(Dir)

    SeedsExtractedFileName=os.path.join("./" + DefenseSystem_Name+"_OUTPUT", "Seeds_" + DefenseSystem_Name + ".tsv")
    VicinityFileName=os.path.join("./" + DefenseSystem_Name+"_OUTPUT", "Vicinity_"+DefenseSystem_Name+".tsv")
    VicinityIDsFileName=os.path.join("./" + DefenseSystem_Name+"_OUTPUT", "VicinityIDs_"+DefenseSystem_Name+".lst")
    VicinityFASTAFileName=os.path.join("./" + DefenseSystem_Name+"_OUTPUT", "Vicinity_"+DefenseSystem_Name+".faa")
    VicinityClustersFileName=os.path.join("./" + DefenseSystem_Name+"_OUTPUT", "VicinityPermissiveClustsLinear_" + DefenseSystem_Name + ".tsv")
    ProfilesFolder=os.path.join("./" + DefenseSystem_Name+"_OUTPUT", "CLUSTERS_"+DefenseSystem_Name + "/")
    SortedBLASTHitsFolder=os.path.join("./" + DefenseSystem_Name+"_OUTPUT", "CLUSTERS_"+DefenseSystem_Name, "Sorted/")
    ICITYFileName=os.path.join("./" + DefenseSystem_Name+"_OUTPUT", "Relevance_"+DefenseSystem_Name+".tsv")
    RelevanceCategoryName="RelevanceCategory_" + DefenseSystem_Name + ".csv"

    # Step 5
    print("Step 5: Extracting seeds")
    HashSeedExtract(PTYFile, SeedPath, SeedsExtractedFileName)

    # Step 6
    FindNeighborhood(PTYFile, SeedsExtractedFileName, VicinityFileName, NeighborhoodVicinitySize)

    # Step 7
    CollectingProteinIDs(VicinityFileName, VicinityIDsFileName)

    # Step 8
    FetchingProteinSequences(PathToDatabase,VicinityIDsFileName,VicinityFASTAFileName)

    # Step 9
    ClusteringProteinSeqiences(VicinityFASTAFileName, PermissiveClusteringThreshold, VicinityClustersFileName)

    # Step 10
    print("Step 10: Making profiles")
    MakeProfiles(VicinityClustersFileName, ProfilesFolder, PathToDatabase)

    # Step 11
    print("Step 11: Running PSI-BLAST for profiles")
    RunPSIBLAST(ProfilesFolder, PathToDatabase, ThreadNum)

    # Step 12
    print("Step 12: Sorting blast hits")
    SortBLASTHitsInMemory(ClustersHitsFolder=ProfilesFolder,
                            SortedHitsFolder=SortedBLASTHitsFolder,
                            PTYFileName=PTYFile,
                            VicinityProteinIDsFileName=VicinityIDsFileName,
                            SeedsFileName=SeedsExtractedFileName,
                            VicinityFileName=VicinityFileName,
                            MIN_OVERLAP_FILTER=SortingOverlapThreshold,
                            MIN_COVERAGE_FILTER=SortingCoverageThresold)

    # Step 13
    print("Step 13: Calculating ICITY metric")
    CalculatingICITYMetric(SortedBLASTHitsFolder=SortedBLASTHitsFolder,
                            PathToDatabase=PathToDatabase,
                            VicinityClustersFileName=VicinityClustersFileName,
                            ICITYFileName=ICITYFileName,
                            ThreadNum=ThreadNum)
