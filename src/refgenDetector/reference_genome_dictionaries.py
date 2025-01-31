#!/usr/bin/env python

""" reference_genome_dictionaries.py: repository with the unique contigs for each reference genome to be inferred"""

__author__ = "Mireia Marin Ginestar"
__version__ = "2.0"
__maintainer__ = "Mireia Marin Ginestar"
__email__ = "mireia.marin@crg.eu"
__status__ = "Developement"


# MAJOR RELEASES
hg16={
    "chr1": 246127941,
    "chr2": 243615958,
    "chr3": 199344050,
    "chr4": 191731959,
    "chr5": 181034922,
    "chr6": 170914576,
    "chr7": 158545518,
    "chr8": 146308819,
    "chr9": 136372045,
    "chrX": 153692391,
    "chrY": 50286555,
    "chr10": 135037215,
    "chr11": 134482954,
    "chr12": 132078379,
    "chr13": 113042980,
    "chr14": 105311216,
    "chr15": 100256656,
    "chr16": 90041932,
    "chr17": 81860266,
    "chr18": 76115139,
    "chr19": 63811651,
    "chr20": 63741868,
    "chr21": 46976097,
    "chr22": 49396972
}

hg17={
    "chr10": 135413628,
    "chr11": 134452384,
    "chr12": 132449811,
    "chr13": 114142980,
    "chr14": 106368585,
    "chr15": 100338915,
    "chr16": 88827254,
    "chr17": 78774742,
    "chr18": 76117153,
    "chr19": 63811651,
    "chr1": 245522847,
    "chr20": 62435964,
    "chr21": 46944323,
    "chr22": 49554710,
    "chr2": 243018229,
    "chr3": 199505740,
    "chr4": 191411218,
    "chr5": 180857866,
    "chr6": 170975699,
    "chr7": 158628139,
    "chr8": 146274826,
    "chr9": 138429268,
    "chrX": 154824264,
    "chrY": 57701691
}

hg18={
    "chr1": 247249719,
    "chr10": 135374737,
    "chr11": 134452384,
    "chr12": 132349534,
    "chr13": 114142980,
    "chr14": 106368585,
    "chr15": 100338915,
    "chr16": 88827254,
    "chr17": 78774742,
    "chr18": 76117153,
    "chr19": 63811651,
    "chr2": 242951149,
    "chr20": 62435964,
    "chr21": 46944323,
    "chr22": 49691432,
    "chr3": 199501827,
    "chr4": 191273063,
    "chr5": 180857866,
    "chr6": 170899992,
    "chr7": 158821424,
    "chr8": 146274826,
    "chr9": 140273252,
    "chrX": 154913754,
    "chrY": 57772954
}
GRCh37={
    "chr1": 249250621,
    "chr2": 243199373,
    "chr3": 198022430,
    "chr4": 191154276,
    "chr5": 180915260,
    "chr6": 171115067,
    "chr7": 159138663,
    "chr8": 146364022,
    "chr9": 141213431,
    "chr10": 135534747,
    "chr11": 135006516,
    "chr12": 133851895,
    "chr13": 115169878,
    "chr14": 107349540,
    "chr15": 102531392,
    "chr16": 90354753,
    "chr17": 81195210,
    "chr18": 78077248,
    "chr19": 59128983,
    "chr20": 63025520,
    "chr21": 48129895,
    "chr22": 51304566,
    "chrX": 155270560,
    "chrY": 59373566
}
GRCh38={
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415
}
T2T={
    "1": 248387328,
    "2": 242696752,
    "3": 201105948,
    "4": 193574945,
    "5": 182045439,
    "6": 172126628,
    "7": 160567428,
    "8": 146259331,
    "9": 150617247,
    "10": 134758134,
    "11": 135127769,
    "12": 133324548,
    "13": 113566686,
    "14": 101161492,
    "15": 99753195,
    "16": 96330374,
    "17": 84276897,
    "18": 80542538,
    "19": 61707364,
    "20": 66210255,
    "21": 45090682,
    "22": 51324926,
    "X": 154259566,
    "Y": 62460029
}

# GRCh37 FLAVORS
hs37d5={"hs37d5":35477943, "NC_007605":171823}
b37={"NC_007605":171823} # this contig is also present in hs37d5, but b37 doesn't include hs37d5 contig
hg19={"ChrM":16571} # chrM it is shared with hg16, hg17 and hg18

# GRCh38 FLAVORS
verily_difGRCh38 = {"chrUn_KN707606v1_decoy": 2200, "chrUn_KN707607v1_decoy": 3033, "chrUn_KN707608v1_decoy": 3112, "chrUn_KN707609v1_decoy": 1642, "chrUn_KN707610v1_decoy": 1393, "chrUn_KN707611v1_decoy": 1103, "chrUn_KN707612v1_decoy": 1039, "chrUn_KN707613v1_decoy": 1619, "chrUn_KN707614v1_decoy": 3122, "chrUn_KN707615v1_decoy": 1934, "chrUn_KN707616v1_decoy": 3111, "chrUn_KN707617v1_decoy": 2545, "chrUn_KN707618v1_decoy": 2295, "chrUn_KN707619v1_decoy": 1551, "chrUn_KN707620v1_decoy": 2046, "chrUn_KN707621v1_decoy": 1222, "chrUn_KN707622v1_decoy": 1535, "chrUn_KN707623v1_decoy": 3784, "chrUn_KN707624v1_decoy": 1329, "chrUn_KN707625v1_decoy": 1238, "chrUn_KN707626v1_decoy": 5623, "chrUn_KN707627v1_decoy": 5821, "chrUn_KN707628v1_decoy": 2960, "chrUn_KN707629v1_decoy": 1848, "chrUn_KN707630v1_decoy": 2315, "chrUn_KN707631v1_decoy": 1945, "chrUn_KN707632v1_decoy": 1424, "chrUn_KN707633v1_decoy": 1274, "chrUn_KN707634v1_decoy": 1007, "chrUn_KN707635v1_decoy": 1414, "chrUn_KN707636v1_decoy": 1725, "chrUn_KN707637v1_decoy": 5354, "chrUn_KN707638v1_decoy": 2189, "chrUn_KN707639v1_decoy": 1294, "chrUn_KN707640v1_decoy": 1831, "chrUn_KN707641v1_decoy": 1647, "chrUn_KN707642v1_decoy": 2943, "chrUn_KN707643v1_decoy": 2857, "chrUn_KN707644v1_decoy": 1030, "chrUn_KN707645v1_decoy": 1070, "chrUn_KN707646v1_decoy": 1735, "chrUn_KN707647v1_decoy": 1982, "chrUn_KN707648v1_decoy": 1564, "chrUn_KN707649v1_decoy": 1775, "chrUn_KN707650v1_decoy": 1540, "chrUn_KN707651v1_decoy": 2013, "chrUn_KN707652v1_decoy": 1176, "chrUn_KN707653v1_decoy": 1890, "chrUn_KN707654v1_decoy": 3644, "chrUn_KN707655v1_decoy": 2785, "chrUn_KN707656v1_decoy": 1017, "chrUn_KN707657v1_decoy": 1068, "chrUn_KN707658v1_decoy": 1007, "chrUn_KN707659v1_decoy": 2605, "chrUn_KN707660v1_decoy": 8410, "chrUn_KN707661v1_decoy": 5534, "chrUn_KN707662v1_decoy": 2173, "chrUn_KN707663v1_decoy": 1065, "chrUn_KN707664v1_decoy": 8683, "chrUn_KN707665v1_decoy": 2670, "chrUn_KN707666v1_decoy": 2420, "chrUn_KN707667v1_decoy": 2189, "chrUn_KN707668v1_decoy": 2093, "chrUn_KN707669v1_decoy": 1184, "chrUn_KN707670v1_decoy": 1205, "chrUn_KN707671v1_decoy": 2786, "chrUn_KN707672v1_decoy": 2794, "chrUn_KN707673v1_decoy": 19544, "chrUn_KN707674v1_decoy": 2848, "chrUn_KN707675v1_decoy": 10556, "chrUn_KN707676v1_decoy": 9066, "chrUn_KN707677v1_decoy": 7267, "chrUn_KN707678v1_decoy": 2462, "chrUn_KN707680v1_decoy": 1297, "chrUn_KN707681v1_decoy": 4379, "chrUn_KN707682v1_decoy": 4208, "chrUn_KN707683v1_decoy": 4068, "chrUn_KN707684v1_decoy": 2940, "chrUn_KN707685v1_decoy": 3938, "chrUn_KN707686v1_decoy": 2072, "chrUn_KN707688v1_decoy": 4248, "chrUn_KN707689v1_decoy": 5823, "chrUn_KN707690v1_decoy": 3715, "chrUn_KN707691v1_decoy": 4885, "chrUn_KN707692v1_decoy": 4813, "chrUn_KN707693v1_decoy": 2899, "chrUn_KN707694v1_decoy": 1228, "chrUn_KN707695v1_decoy": 3119, "chrUn_KN707696v1_decoy": 3828, "chrUn_KN707697v1_decoy": 1186, "chrUn_KN707698v1_decoy": 1908, "chrUn_KN707699v1_decoy": 2795, "chrUn_KN707700v1_decoy": 3703, "chrUn_KN707701v1_decoy": 6722, "chrUn_KN707702v1_decoy": 6466, "chrUn_KN707703v1_decoy": 2235, "chrUn_KN707704v1_decoy": 2871, "chrUn_KN707705v1_decoy": 4632, "chrUn_KN707706v1_decoy": 4225, "chrUn_KN707707v1_decoy": 4339, "chrUn_KN707708v1_decoy": 2305, "chrUn_KN707709v1_decoy": 3273, "chrUn_KN707710v1_decoy": 5701, "chrUn_KN707711v1_decoy": 4154, "chrUn_KN707712v1_decoy": 1243, "chrUn_KN707714v1_decoy": 2922, "chrUn_KN707715v1_decoy": 3044, "chrUn_KN707716v1_decoy": 2888, "chrUn_KN707717v1_decoy": 1742, "chrUn_KN707718v1_decoy": 4969, "chrUn_KN707719v1_decoy": 3270, "chrUn_KN707720v1_decoy": 6028, "chrUn_KN707721v1_decoy": 1105, "chrUn_KN707722v1_decoy": 2884, "chrUn_KN707723v1_decoy": 1124, "chrUn_KN707724v1_decoy": 1454, "chrUn_KN707725v1_decoy": 2565, "chrUn_KN707726v1_decoy": 2149, "chrUn_KN707727v1_decoy": 2630, "chrUn_KN707728v1_decoy": 14625, "chrUn_KN707729v1_decoy": 7431, "chrUn_KN707730v1_decoy": 5776, "chrUn_KN707731v1_decoy": 4820, "chrUn_KN707732v1_decoy": 1227, "chrUn_KN707733v1_decoy": 7503, "chrUn_KN707734v1_decoy": 9652, "chrUn_KN707735v1_decoy": 1091, "chrUn_KN707736v1_decoy": 2467, "chrUn_KN707737v1_decoy": 1270, "chrUn_KN707738v1_decoy": 4365, "chrUn_KN707739v1_decoy": 4284, "chrUn_KN707740v1_decoy": 10282, "chrUn_KN707741v1_decoy": 5601, "chrUn_KN707742v1_decoy": 4758, "chrUn_KN707743v1_decoy": 1624, "chrUn_KN707744v1_decoy": 4024, "chrUn_KN707745v1_decoy": 1276, "chrUn_KN707746v1_decoy": 5083, "chrUn_KN707747v1_decoy": 2075, "chrUn_KN707748v1_decoy": 3553, "chrUn_KN707749v1_decoy": 7010, "chrUn_KN707750v1_decoy": 4718, "chrUn_KN707751v1_decoy": 3546, "chrUn_KN707752v1_decoy": 2873, "chrUn_KN707753v1_decoy": 2144, "chrUn_KN707754v1_decoy": 2243, "chrUn_KN707755v1_decoy": 5343, "chrUn_KN707756v1_decoy": 4877, "chrUn_KN707757v1_decoy": 3034, "chrUn_KN707758v1_decoy": 2826, "chrUn_KN707759v1_decoy": 1221, "chrUn_KN707760v1_decoy": 1169, "chrUn_KN707761v1_decoy": 2319, "chrUn_KN707762v1_decoy": 3450, "chrUn_KN707763v1_decoy": 2674, "chrUn_KN707764v1_decoy": 3912, "chrUn_KN707765v1_decoy": 6020, "chrUn_KN707766v1_decoy": 2303, "chrUn_KN707767v1_decoy": 2552, "chrUn_KN707768v1_decoy": 3656, "chrUn_KN707769v1_decoy": 1591, "chrUn_KN707770v1_decoy": 1209, "chrUn_KN707771v1_decoy": 3176, "chrUn_KN707772v1_decoy": 8915, "chrUn_KN707773v1_decoy": 4902, "chrUn_KN707774v1_decoy": 3324, "chrUn_KN707775v1_decoy": 5997, "chrUn_KN707776v1_decoy": 2618, "chrUn_KN707777v1_decoy": 10311, "chrUn_KN707778v1_decoy": 2440, "chrUn_KN707779v1_decoy": 12444, "chrUn_KN707780v1_decoy": 5691, "chrUn_KN707781v1_decoy": 2717, "chrUn_KN707782v1_decoy": 5277, "chrUn_KN707783v1_decoy": 4373, "chrUn_KN707784v1_decoy": 3224, "chrUn_KN707785v1_decoy": 2631, "chrUn_KN707786v1_decoy": 5385, "chrUn_KN707787v1_decoy": 3678, "chrUn_KN707788v1_decoy": 1412, "chrUn_KN707789v1_decoy": 1443, "chrUn_KN707790v1_decoy": 1098, "chrUn_KN707791v1_decoy": 3240, "chrUn_KN707792v1_decoy": 1915, "chrUn_KN707793v1_decoy": 4667, "chrUn_KN707794v1_decoy": 7219, "chrUn_KN707795v1_decoy": 3277, "chrUn_KN707796v1_decoy": 3473, "chrUn_KN707797v1_decoy": 4243, "chrUn_KN707798v1_decoy": 17599, "chrUn_KN707799v1_decoy": 5095, "chrUn_KN707800v1_decoy": 2237, "chrUn_KN707801v1_decoy": 2901, "chrUn_KN707802v1_decoy": 2666, "chrUn_KN707803v1_decoy": 5336, "chrUn_KN707804v1_decoy": 4383, "chrUn_KN707805v1_decoy": 5446, "chrUn_KN707806v1_decoy": 6252, "chrUn_KN707807v1_decoy": 4616, "chrUn_KN707808v1_decoy": 3021, "chrUn_KN707809v1_decoy": 3667, "chrUn_KN707810v1_decoy": 4563, "chrUn_KN707811v1_decoy": 1120, "chrUn_KN707812v1_decoy": 3845, "chrUn_KN707813v1_decoy": 2272, "chrUn_KN707814v1_decoy": 4764, "chrUn_KN707815v1_decoy": 5410, "chrUn_KN707816v1_decoy": 7150, "chrUn_KN707817v1_decoy": 1762, "chrUn_KN707818v1_decoy": 1207, "chrUn_KN707819v1_decoy": 1331, "chrUn_KN707820v1_decoy": 8307, "chrUn_KN707822v1_decoy": 2575, "chrUn_KN707823v1_decoy": 3970, "chrUn_KN707824v1_decoy": 1352, "chrUn_KN707825v1_decoy": 3040, "chrUn_KN707826v1_decoy": 2070, "chrUn_KN707827v1_decoy": 2913, "chrUn_KN707828v1_decoy": 2389, "chrUn_KN707829v1_decoy": 1835, "chrUn_KN707830v1_decoy": 4807, "chrUn_KN707831v1_decoy": 2201, "chrUn_KN707832v1_decoy": 1265, "chrUn_KN707833v1_decoy": 1961, "chrUn_KN707834v1_decoy": 1064, "chrUn_KN707835v1_decoy": 1932, "chrUn_KN707836v1_decoy": 3213, "chrUn_KN707837v1_decoy": 1178, "chrUn_KN707838v1_decoy": 2926, "chrUn_KN707839v1_decoy": 1038, "chrUn_KN707840v1_decoy": 3298, "chrUn_KN707841v1_decoy": 8992, "chrUn_KN707842v1_decoy": 6698, "chrUn_KN707843v1_decoy": 4880, "chrUn_KN707844v1_decoy": 1766, "chrUn_KN707845v1_decoy": 3532, "chrUn_KN707846v1_decoy": 2297, "chrUn_KN707847v1_decoy": 1234, "chrUn_KN707848v1_decoy": 1205, "chrUn_KN707849v1_decoy": 2790, "chrUn_KN707850v1_decoy": 2006, "chrUn_KN707851v1_decoy": 4593, "chrUn_KN707852v1_decoy": 1579, "chrUn_KN707853v1_decoy": 9597, "chrUn_KN707854v1_decoy": 10451, "chrUn_KN707855v1_decoy": 3219, "chrUn_KN707856v1_decoy": 2300, "chrUn_KN707857v1_decoy": 5985, "chrUn_KN707858v1_decoy": 2959, "chrUn_KN707859v1_decoy": 1340, "chrUn_KN707860v1_decoy": 3148, "chrUn_KN707861v1_decoy": 2242, "chrUn_KN707862v1_decoy": 16513, "chrUn_KN707863v1_decoy": 7821, "chrUn_KN707864v1_decoy": 2159, "chrUn_KN707865v1_decoy": 2114, "chrUn_KN707866v1_decoy": 4109, "chrUn_KN707867v1_decoy": 1544, "chrUn_KN707868v1_decoy": 1005, "chrUn_KN707869v1_decoy": 8632, "chrUn_KN707870v1_decoy": 1012, "chrUn_KN707871v1_decoy": 4728, "chrUn_KN707873v1_decoy": 7591, "chrUn_KN707874v1_decoy": 5202, "chrUn_KN707875v1_decoy": 4241, "chrUn_KN707876v1_decoy": 4131, "chrUn_KN707877v1_decoy": 2272, "chrUn_KN707878v1_decoy": 2085, "chrUn_KN707879v1_decoy": 4346, "chrUn_KN707880v1_decoy": 1208, "chrUn_KN707881v1_decoy": 4543, "chrUn_KN707882v1_decoy": 2772, "chrUn_KN707883v1_decoy": 2490, "chrUn_KN707884v1_decoy": 4568, "chrUn_KN707885v1_decoy": 1776, "chrUn_KN707887v1_decoy": 3534, "chrUn_KN707888v1_decoy": 2424, "chrUn_KN707889v1_decoy": 1747, "chrUn_KN707890v1_decoy": 1088, "chrUn_KN707892v1_decoy": 2530, "chrUn_KN707893v1_decoy": 8049, "chrUn_KN707894v1_decoy": 1366, "chrUn_KN707895v1_decoy": 4284, "chrUn_KN707896v1_decoy": 33125, "chrUn_KN707897v1_decoy": 2137, "chrUn_KN707898v1_decoy": 3840, "chrUn_KN707899v1_decoy": 3087, "chrUn_KN707900v1_decoy": 2041, "chrUn_KN707901v1_decoy": 3344, "chrUn_KN707902v1_decoy": 2921, "chrUn_KN707903v1_decoy": 6581, "chrUn_KN707904v1_decoy": 3968, "chrUn_KN707905v1_decoy": 2339, "chrUn_KN707906v1_decoy": 1243, "chrUn_KN707907v1_decoy": 7776, "chrUn_KN707908v1_decoy": 19837, "chrUn_KN707909v1_decoy": 1737, "chrUn_KN707910v1_decoy": 1098, "chrUn_KN707911v1_decoy": 1893, "chrUn_KN707912v1_decoy": 1281, "chrUn_KN707913v1_decoy": 1527, "chrUn_KN707914v1_decoy": 2055, "chrUn_KN707915v1_decoy": 2527, "chrUn_KN707916v1_decoy": 3275, "chrUn_KN707917v1_decoy": 1265, "chrUn_KN707918v1_decoy": 2623, "chrUn_KN707919v1_decoy": 4850, "chrUn_KN707920v1_decoy": 3584, "chrUn_KN707921v1_decoy": 2561, "chrUn_KN707923v1_decoy": 1409, "chrUn_KN707924v1_decoy": 4596, "chrUn_KN707925v1_decoy": 11555, "chrUn_KN707926v1_decoy": 1266, "chrUn_KN707927v1_decoy": 1079, "chrUn_KN707928v1_decoy": 1087, "chrUn_KN707929v1_decoy": 1226, "chrUn_KN707930v1_decoy": 1131, "chrUn_KN707931v1_decoy": 1199, "chrUn_KN707932v1_decoy": 1084, "chrUn_KN707933v1_decoy": 2038, "chrUn_KN707934v1_decoy": 1070, "chrUn_KN707935v1_decoy": 1312, "chrUn_KN707936v1_decoy": 4031, "chrUn_KN707937v1_decoy": 7445, "chrUn_KN707938v1_decoy": 1770, "chrUn_KN707939v1_decoy": 5600, "chrUn_KN707940v1_decoy": 1882, "chrUn_KN707941v1_decoy": 1170, "chrUn_KN707943v1_decoy": 5325, "chrUn_KN707945v1_decoy": 1072, "chrUn_KN707946v1_decoy": 2463, "chrUn_KN707947v1_decoy": 1010, "chrUn_KN707948v1_decoy": 1432, "chrUn_KN707949v1_decoy": 1162, "chrUn_KN707950v1_decoy": 1095, "chrUn_KN707951v1_decoy": 1118, "chrUn_KN707952v1_decoy": 1383, "chrUn_KN707953v1_decoy": 2289, "chrUn_KN707954v1_decoy": 1648, "chrUn_KN707955v1_decoy": 2203, "chrUn_KN707956v1_decoy": 3270, "chrUn_KN707957v1_decoy": 11499, "chrUn_KN707958v1_decoy": 2474, "chrUn_KN707959v1_decoy": 2294, "chrUn_KN707960v1_decoy": 1238, "chrUn_KN707961v1_decoy": 3410, "chrUn_KN707962v1_decoy": 1523, "chrUn_KN707963v1_decoy": 62955, "chrUn_KN707964v1_decoy": 6282, "chrUn_KN707965v1_decoy": 3836, "chrUn_KN707966v1_decoy": 6486, "chrUn_KN707967v1_decoy": 15368, "chrUn_KN707968v1_decoy": 9572, "chrUn_KN707969v1_decoy": 6413, "chrUn_KN707970v1_decoy": 4104, "chrUn_KN707971v1_decoy": 12943, "chrUn_KN707972v1_decoy": 4650, "chrUn_KN707973v1_decoy": 3080, "chrUn_KN707974v1_decoy": 3134, "chrUn_KN707975v1_decoy": 6211, "chrUn_KN707976v1_decoy": 1126, "chrUn_KN707977v1_decoy": 1101, "chrUn_KN707978v1_decoy": 1101, "chrUn_KN707979v1_decoy": 2648, "chrUn_KN707980v1_decoy": 2973, "chrUn_KN707981v1_decoy": 2520, "chrUn_KN707983v1_decoy": 2606, "chrUn_KN707984v1_decoy": 2205, "chrUn_KN707985v1_decoy": 2929, "chrUn_KN707986v1_decoy": 3869, "chrUn_KN707987v1_decoy": 1117, "chrUn_KN707988v1_decoy": 2960, "chrUn_KN707989v1_decoy": 1009, "chrUn_KN707990v1_decoy": 4048, "chrUn_KN707991v1_decoy": 2193, "chrUn_KN707992v1_decoy": 1830, "chrUn_JTFH01000001v1_decoy": 25139, "chrUn_JTFH01000002v1_decoy": 18532, "chrUn_JTFH01000003v1_decoy": 15240, "chrUn_JTFH01000004v1_decoy": 13739, "chrUn_JTFH01000005v1_decoy": 11297, "chrUn_JTFH01000006v1_decoy": 10074, "chrUn_JTFH01000007v1_decoy": 9891, "chrUn_JTFH01000008v1_decoy": 9774, "chrUn_JTFH01000009v1_decoy": 9727, "chrUn_JTFH01000010v1_decoy": 9358, "chrUn_JTFH01000011v1_decoy": 8920, "chrUn_JTFH01000012v1_decoy": 8479, "chrUn_JTFH01000013v1_decoy": 8312, "chrUn_JTFH01000014v1_decoy": 8261, "chrUn_JTFH01000015v1_decoy": 8131, "chrUn_JTFH01000016v1_decoy": 8051, "chrUn_JTFH01000017v1_decoy": 7832, "chrUn_JTFH01000018v1_decoy": 7710, "chrUn_JTFH01000019v1_decoy": 7702, "chrUn_JTFH01000020v1_decoy": 7479, "chrUn_JTFH01000021v1_decoy": 7368, "chrUn_JTFH01000022v1_decoy": 7162, "chrUn_JTFH01000023v1_decoy": 7065, "chrUn_JTFH01000024v1_decoy": 7019, "chrUn_JTFH01000025v1_decoy": 6997, "chrUn_JTFH01000026v1_decoy": 6994, "chrUn_JTFH01000027v1_decoy": 6979, "chrUn_JTFH01000028v1_decoy": 6797, "chrUn_JTFH01000029v1_decoy": 6525, "chrUn_JTFH01000030v1_decoy": 6246, "chrUn_JTFH01000031v1_decoy": 5926, "chrUn_JTFH01000032v1_decoy": 5914, "chrUn_JTFH01000033v1_decoy": 5898, "chrUn_JTFH01000034v1_decoy": 5879, "chrUn_JTFH01000035v1_decoy": 5834, "chrUn_JTFH01000036v1_decoy": 5743, "chrUn_JTFH01000037v1_decoy": 5577, "chrUn_JTFH01000038v1_decoy": 5413, "chrUn_JTFH01000039v1_decoy": 5250, "chrUn_JTFH01000040v1_decoy": 5246, "chrUn_JTFH01000041v1_decoy": 5118, "chrUn_JTFH01000042v1_decoy": 5058, "chrUn_JTFH01000043v1_decoy": 4959, "chrUn_JTFH01000044v1_decoy": 4853, "chrUn_JTFH01000045v1_decoy": 4828, "chrUn_JTFH01000046v1_decoy": 4819, "chrUn_JTFH01000047v1_decoy": 4809, "chrUn_JTFH01000048v1_decoy": 4710, "chrUn_JTFH01000049v1_decoy": 4680, "chrUn_JTFH01000050v1_decoy": 4645, "chrUn_JTFH01000051v1_decoy": 4514, "chrUn_JTFH01000052v1_decoy": 4439, "chrUn_JTFH01000054v1_decoy": 4409, "chrUn_JTFH01000055v1_decoy": 4392, "chrUn_JTFH01000056v1_decoy": 4359, "chrUn_JTFH01000057v1_decoy": 4319, "chrUn_JTFH01000058v1_decoy": 4290, "chrUn_JTFH01000059v1_decoy": 4242, "chrUn_JTFH01000060v1_decoy": 4228, "chrUn_JTFH01000061v1_decoy": 4222, "chrUn_JTFH01000062v1_decoy": 4216, "chrUn_JTFH01000063v1_decoy": 4210, "chrUn_JTFH01000064v1_decoy": 4206, "chrUn_JTFH01000065v1_decoy": 4102, "chrUn_JTFH01000066v1_decoy": 4101, "chrUn_JTFH01000067v1_decoy": 4083, "chrUn_JTFH01000068v1_decoy": 3967, "chrUn_JTFH01000069v1_decoy": 3955, "chrUn_JTFH01000070v1_decoy": 3945, "chrUn_JTFH01000071v1_decoy": 3930, "chrUn_JTFH01000072v1_decoy": 3929, "chrUn_JTFH01000073v1_decoy": 3924, "chrUn_JTFH01000074v1_decoy": 3919, "chrUn_JTFH01000075v1_decoy": 3908, "chrUn_JTFH01000076v1_decoy": 3892, "chrUn_JTFH01000077v1_decoy": 3890, "chrUn_JTFH01000078v1_decoy": 3859, "chrUn_JTFH01000079v1_decoy": 3846, "chrUn_JTFH01000080v1_decoy": 3835, "chrUn_JTFH01000081v1_decoy": 3830, "chrUn_JTFH01000082v1_decoy": 3828, "chrUn_JTFH01000083v1_decoy": 3825, "chrUn_JTFH01000084v1_decoy": 3821, "chrUn_JTFH01000085v1_decoy": 3809, "chrUn_JTFH01000086v1_decoy": 3801, "chrUn_JTFH01000087v1_decoy": 3799, "chrUn_JTFH01000088v1_decoy": 3737, "chrUn_JTFH01000089v1_decoy": 3701, "chrUn_JTFH01000090v1_decoy": 3698, "chrUn_JTFH01000091v1_decoy": 3692, "chrUn_JTFH01000092v1_decoy": 3686, "chrUn_JTFH01000093v1_decoy": 3677, "chrUn_JTFH01000094v1_decoy": 3664, "chrUn_JTFH01000095v1_decoy": 3613, "chrUn_JTFH01000096v1_decoy": 3611, "chrUn_JTFH01000097v1_decoy": 3606, "chrUn_JTFH01000098v1_decoy": 3584, "chrUn_JTFH01000099v1_decoy": 3581, "chrUn_JTFH01000100v1_decoy": 3543, "chrUn_JTFH01000101v1_decoy": 3528, "chrUn_JTFH01000102v1_decoy": 3527, "chrUn_JTFH01000103v1_decoy": 3496, "chrUn_JTFH01000104v1_decoy": 3493, "chrUn_JTFH01000105v1_decoy": 3484, "chrUn_JTFH01000106v1_decoy": 3435, "chrUn_JTFH01000107v1_decoy": 3391, "chrUn_JTFH01000108v1_decoy": 3374, "chrUn_JTFH01000109v1_decoy": 3371, "chrUn_JTFH01000110v1_decoy": 3361, "chrUn_JTFH01000111v1_decoy": 3351, "chrUn_JTFH01000112v1_decoy": 3340, "chrUn_JTFH01000113v1_decoy": 3320, "chrUn_JTFH01000114v1_decoy": 3282, "chrUn_JTFH01000115v1_decoy": 3278, "chrUn_JTFH01000116v1_decoy": 3260, "chrUn_JTFH01000117v1_decoy": 3258, "chrUn_JTFH01000119v1_decoy": 3247, "chrUn_JTFH01000120v1_decoy": 3230, "chrUn_JTFH01000121v1_decoy": 3224, "chrUn_JTFH01000122v1_decoy": 3216, "chrUn_JTFH01000123v1_decoy": 3212, "chrUn_JTFH01000124v1_decoy": 3194, "chrUn_JTFH01000125v1_decoy": 3189, "chrUn_JTFH01000126v1_decoy": 3177, "chrUn_JTFH01000127v1_decoy": 3176, "chrUn_JTFH01000128v1_decoy": 3173, "chrUn_JTFH01000129v1_decoy": 3170, "chrUn_JTFH01000130v1_decoy": 3166, "chrUn_JTFH01000131v1_decoy": 3163, "chrUn_JTFH01000132v1_decoy": 3143, "chrUn_JTFH01000133v1_decoy": 3137, "chrUn_JTFH01000134v1_decoy": 3116, "chrUn_JTFH01000135v1_decoy": 3106, "chrUn_JTFH01000136v1_decoy": 3093, "chrUn_JTFH01000137v1_decoy": 3079, "chrUn_JTFH01000138v1_decoy": 3053, "chrUn_JTFH01000139v1_decoy": 3051, "chrUn_JTFH01000140v1_decoy": 3015, "chrUn_JTFH01000141v1_decoy": 3012, "chrUn_JTFH01000142v1_decoy": 3009, "chrUn_JTFH01000143v1_decoy": 2997, "chrUn_JTFH01000144v1_decoy": 2997, "chrUn_JTFH01000146v1_decoy": 2979, "chrUn_JTFH01000147v1_decoy": 2967, "chrUn_JTFH01000148v1_decoy": 2967, "chrUn_JTFH01000149v1_decoy": 2966, "chrUn_JTFH01000150v1_decoy": 2954, "chrUn_JTFH01000151v1_decoy": 2952, "chrUn_JTFH01000152v1_decoy": 2934, "chrUn_JTFH01000153v1_decoy": 2918, "chrUn_JTFH01000154v1_decoy": 2895, "chrUn_JTFH01000155v1_decoy": 2887, "chrUn_JTFH01000156v1_decoy": 2879, "chrUn_JTFH01000157v1_decoy": 2878, "chrUn_JTFH01000158v1_decoy": 2872, "chrUn_JTFH01000159v1_decoy": 2868, "chrUn_JTFH01000160v1_decoy": 2866, "chrUn_JTFH01000161v1_decoy": 2865, "chrUn_JTFH01000162v1_decoy": 2864, "chrUn_JTFH01000163v1_decoy": 2859, "chrUn_JTFH01000164v1_decoy": 2854, "chrUn_JTFH01000165v1_decoy": 2830, "chrUn_JTFH01000166v1_decoy": 2828, "chrUn_JTFH01000167v1_decoy": 2824, "chrUn_JTFH01000168v1_decoy": 2819, "chrUn_JTFH01000169v1_decoy": 2813, "chrUn_JTFH01000170v1_decoy": 2809, "chrUn_JTFH01000171v1_decoy": 2802, "chrUn_JTFH01000172v1_decoy": 2791, "chrUn_JTFH01000173v1_decoy": 2783, "chrUn_JTFH01000174v1_decoy": 2778, "chrUn_JTFH01000175v1_decoy": 2777, "chrUn_JTFH01000176v1_decoy": 2770, "chrUn_JTFH01000177v1_decoy": 2769, "chrUn_JTFH01000178v1_decoy": 2767, "chrUn_JTFH01000179v1_decoy": 2763, "chrUn_JTFH01000180v1_decoy": 2745, "chrUn_JTFH01000181v1_decoy": 2742, "chrUn_JTFH01000182v1_decoy": 2736, "chrUn_JTFH01000183v1_decoy": 2729, "chrUn_JTFH01000184v1_decoy": 2726, "chrUn_JTFH01000185v1_decoy": 2719, "chrUn_JTFH01000186v1_decoy": 2715, "chrUn_JTFH01000187v1_decoy": 2708, "chrUn_JTFH01000188v1_decoy": 2704, "chrUn_JTFH01000189v1_decoy": 2692, "chrUn_JTFH01000190v1_decoy": 2691, "chrUn_JTFH01000191v1_decoy": 2690, "chrUn_JTFH01000192v1_decoy": 2687, "chrUn_JTFH01000193v1_decoy": 2677, "chrUn_JTFH01000194v1_decoy": 2668, "chrUn_JTFH01000195v1_decoy": 2668, "chrUn_JTFH01000196v1_decoy": 2663, "chrUn_JTFH01000197v1_decoy": 2655, "chrUn_JTFH01000198v1_decoy": 2644, "chrUn_JTFH01000199v1_decoy": 2642, "chrUn_JTFH01000200v1_decoy": 2632, "chrUn_JTFH01000201v1_decoy": 2632, "chrUn_JTFH01000202v1_decoy": 2628, "chrUn_JTFH01000203v1_decoy": 2623, "chrUn_JTFH01000204v1_decoy": 2622, "chrUn_JTFH01000205v1_decoy": 2619, "chrUn_JTFH01000206v1_decoy": 2605, "chrUn_JTFH01000207v1_decoy": 2603, "chrUn_JTFH01000208v1_decoy": 2601, "chrUn_JTFH01000209v1_decoy": 2598, "chrUn_JTFH01000210v1_decoy": 2597, "chrUn_JTFH01000211v1_decoy": 2596, "chrUn_JTFH01000212v1_decoy": 2594, "chrUn_JTFH01000213v1_decoy": 2586, "chrUn_JTFH01000214v1_decoy": 2585, "chrUn_JTFH01000215v1_decoy": 2583, "chrUn_JTFH01000216v1_decoy": 2578, "chrUn_JTFH01000217v1_decoy": 2569, "chrUn_JTFH01000218v1_decoy": 2569, "chrUn_JTFH01000219v1_decoy": 2551, "chrUn_JTFH01000220v1_decoy": 2548, "chrUn_JTFH01000221v1_decoy": 2548, "chrUn_JTFH01000222v1_decoy": 2546, "chrUn_JTFH01000223v1_decoy": 2545, "chrUn_JTFH01000224v1_decoy": 2534, "chrUn_JTFH01000225v1_decoy": 2533, "chrUn_JTFH01000226v1_decoy": 2522, "chrUn_JTFH01000227v1_decoy": 2522, "chrUn_JTFH01000228v1_decoy": 2515, "chrUn_JTFH01000229v1_decoy": 2513, "chrUn_JTFH01000230v1_decoy": 2507, "chrUn_JTFH01000231v1_decoy": 2504, "chrUn_JTFH01000232v1_decoy": 2497, "chrUn_JTFH01000233v1_decoy": 2471, "chrUn_JTFH01000234v1_decoy": 2465, "chrUn_JTFH01000235v1_decoy": 2464, "chrUn_JTFH01000236v1_decoy": 2459, "chrUn_JTFH01000237v1_decoy": 2457, "chrUn_JTFH01000238v1_decoy": 2450, "chrUn_JTFH01000239v1_decoy": 2435, "chrUn_JTFH01000240v1_decoy": 2434, "chrUn_JTFH01000241v1_decoy": 2432, "chrUn_JTFH01000242v1_decoy": 2427, "chrUn_JTFH01000243v1_decoy": 2421, "chrUn_JTFH01000244v1_decoy": 2420, "chrUn_JTFH01000245v1_decoy": 2414, "chrUn_JTFH01000246v1_decoy": 2404, "chrUn_JTFH01000247v1_decoy": 2403, "chrUn_JTFH01000248v1_decoy": 2402, "chrUn_JTFH01000249v1_decoy": 2397, "chrUn_JTFH01000250v1_decoy": 2395, "chrUn_JTFH01000251v1_decoy": 2394, "chrUn_JTFH01000252v1_decoy": 2388, "chrUn_JTFH01000253v1_decoy": 2382, "chrUn_JTFH01000254v1_decoy": 2381, "chrUn_JTFH01000255v1_decoy": 2380, "chrUn_JTFH01000256v1_decoy": 2368, "chrUn_JTFH01000257v1_decoy": 2364, "chrUn_JTFH01000258v1_decoy": 2363, "chrUn_JTFH01000259v1_decoy": 2348, "chrUn_JTFH01000260v1_decoy": 2339, "chrUn_JTFH01000261v1_decoy": 2335, "chrUn_JTFH01000262v1_decoy": 2332, "chrUn_JTFH01000263v1_decoy": 2331, "chrUn_JTFH01000264v1_decoy": 2330, "chrUn_JTFH01000265v1_decoy": 2323, "chrUn_JTFH01000266v1_decoy": 2319, "chrUn_JTFH01000267v1_decoy": 2314, "chrUn_JTFH01000268v1_decoy": 2308, "chrUn_JTFH01000269v1_decoy": 2306, "chrUn_JTFH01000270v1_decoy": 2296, "chrUn_JTFH01000271v1_decoy": 2287, "chrUn_JTFH01000272v1_decoy": 2279, "chrUn_JTFH01000274v1_decoy": 2273, "chrUn_JTFH01000275v1_decoy": 2262, "chrUn_JTFH01000276v1_decoy": 2254, "chrUn_JTFH01000277v1_decoy": 2252, "chrUn_JTFH01000278v1_decoy": 2245, "chrUn_JTFH01000279v1_decoy": 2239, "chrUn_JTFH01000280v1_decoy": 2223, "chrUn_JTFH01000281v1_decoy": 2220, "chrUn_JTFH01000282v1_decoy": 2218, "chrUn_JTFH01000283v1_decoy": 2215, "chrUn_JTFH01000284v1_decoy": 2213, "chrUn_JTFH01000285v1_decoy": 2203, "chrUn_JTFH01000286v1_decoy": 2200, "chrUn_JTFH01000287v1_decoy": 2197, "chrUn_JTFH01000288v1_decoy": 2194, "chrUn_JTFH01000289v1_decoy": 2183, "chrUn_JTFH01000290v1_decoy": 2179, "chrUn_JTFH01000291v1_decoy": 2177, "chrUn_JTFH01000292v1_decoy": 2177, "chrUn_JTFH01000293v1_decoy": 2177, "chrUn_JTFH01000295v1_decoy": 2160, "chrUn_JTFH01000296v1_decoy": 2155, "chrUn_JTFH01000297v1_decoy": 2144, "chrUn_JTFH01000298v1_decoy": 2143, "chrUn_JTFH01000299v1_decoy": 2136, "chrUn_JTFH01000300v1_decoy": 2134, "chrUn_JTFH01000301v1_decoy": 2129, "chrUn_JTFH01000302v1_decoy": 2128, "chrUn_JTFH01000303v1_decoy": 2125, "chrUn_JTFH01000304v1_decoy": 2125, "chrUn_JTFH01000305v1_decoy": 2122, "chrUn_JTFH01000306v1_decoy": 2111, "chrUn_JTFH01000307v1_decoy": 2106, "chrUn_JTFH01000308v1_decoy": 2094, "chrUn_JTFH01000309v1_decoy": 2093, "chrUn_JTFH01000310v1_decoy": 2088, "chrUn_JTFH01000311v1_decoy": 2086, "chrUn_JTFH01000312v1_decoy": 2086, "chrUn_JTFH01000313v1_decoy": 2084, "chrUn_JTFH01000314v1_decoy": 2080, "chrUn_JTFH01000315v1_decoy": 2079, "chrUn_JTFH01000316v1_decoy": 2076, "chrUn_JTFH01000317v1_decoy": 2071, "chrUn_JTFH01000318v1_decoy": 2066, "chrUn_JTFH01000319v1_decoy": 2061, "chrUn_JTFH01000320v1_decoy": 2055, "chrUn_JTFH01000321v1_decoy": 2053, "chrUn_JTFH01000322v1_decoy": 2040, "chrUn_JTFH01000323v1_decoy": 2036, "chrUn_JTFH01000324v1_decoy": 2035, "chrUn_JTFH01000325v1_decoy": 2034, "chrUn_JTFH01000326v1_decoy": 2032, "chrUn_JTFH01000327v1_decoy": 2029, "chrUn_JTFH01000328v1_decoy": 2025, "chrUn_JTFH01000329v1_decoy": 2021, "chrUn_JTFH01000330v1_decoy": 2018, "chrUn_JTFH01000331v1_decoy": 2015, "chrUn_JTFH01000332v1_decoy": 2009, "chrUn_JTFH01000333v1_decoy": 2007, "chrUn_JTFH01000334v1_decoy": 2005, "chrUn_JTFH01000335v1_decoy": 2003, "chrUn_JTFH01000336v1_decoy": 2001, "chrUn_JTFH01000337v1_decoy": 2001, "chrUn_JTFH01000338v1_decoy": 2000, "chrUn_JTFH01000339v1_decoy": 1996, "chrUn_JTFH01000340v1_decoy": 1992, "chrUn_JTFH01000341v1_decoy": 1985, "chrUn_JTFH01000342v1_decoy": 1981, "chrUn_JTFH01000343v1_decoy": 1977, "chrUn_JTFH01000344v1_decoy": 1971, "chrUn_JTFH01000345v1_decoy": 1968, "chrUn_JTFH01000346v1_decoy": 1962, "chrUn_JTFH01000347v1_decoy": 1961, "chrUn_JTFH01000348v1_decoy": 1960, "chrUn_JTFH01000349v1_decoy": 1960, "chrUn_JTFH01000350v1_decoy": 1954, "chrUn_JTFH01000351v1_decoy": 1952, "chrUn_JTFH01000352v1_decoy": 1947, "chrUn_JTFH01000353v1_decoy": 1944, "chrUn_JTFH01000354v1_decoy": 1943, "chrUn_JTFH01000355v1_decoy": 1941, "chrUn_JTFH01000356v1_decoy": 1937, "chrUn_JTFH01000357v1_decoy": 1934, "chrUn_JTFH01000358v1_decoy": 1929, "chrUn_JTFH01000359v1_decoy": 1924, "chrUn_JTFH01000360v1_decoy": 1924, "chrUn_JTFH01000361v1_decoy": 1923, "chrUn_JTFH01000362v1_decoy": 1921, "chrUn_JTFH01000363v1_decoy": 1918, "chrUn_JTFH01000364v1_decoy": 1915, "chrUn_JTFH01000365v1_decoy": 1915, "chrUn_JTFH01000366v1_decoy": 1914, "chrUn_JTFH01000367v1_decoy": 1912, "chrUn_JTFH01000368v1_decoy": 1910, "chrUn_JTFH01000369v1_decoy": 1907, "chrUn_JTFH01000370v1_decoy": 1904, "chrUn_JTFH01000371v1_decoy": 1897, "chrUn_JTFH01000372v1_decoy": 1891, "chrUn_JTFH01000373v1_decoy": 1890, "chrUn_JTFH01000374v1_decoy": 1888, "chrUn_JTFH01000375v1_decoy": 1888, "chrUn_JTFH01000376v1_decoy": 1885, "chrUn_JTFH01000377v1_decoy": 1881, "chrUn_JTFH01000378v1_decoy": 1879, "chrUn_JTFH01000379v1_decoy": 1877, "chrUn_JTFH01000380v1_decoy": 1876, "chrUn_JTFH01000381v1_decoy": 1876, "chrUn_JTFH01000382v1_decoy": 1874, "chrUn_JTFH01000383v1_decoy": 1872, "chrUn_JTFH01000384v1_decoy": 1869, "chrUn_JTFH01000385v1_decoy": 1866, "chrUn_JTFH01000386v1_decoy": 1865, "chrUn_JTFH01000387v1_decoy": 1865, "chrUn_JTFH01000388v1_decoy": 1865, "chrUn_JTFH01000389v1_decoy": 1862, "chrUn_JTFH01000390v1_decoy": 1862, "chrUn_JTFH01000391v1_decoy": 1859, "chrUn_JTFH01000392v1_decoy": 1856, "chrUn_JTFH01000393v1_decoy": 1856, "chrUn_JTFH01000394v1_decoy": 1854, "chrUn_JTFH01000395v1_decoy": 1850, "chrUn_JTFH01000396v1_decoy": 1849, "chrUn_JTFH01000397v1_decoy": 1849, "chrUn_JTFH01000398v1_decoy": 1847, "chrUn_JTFH01000399v1_decoy": 1839, "chrUn_JTFH01000400v1_decoy": 1834, "chrUn_JTFH01000401v1_decoy": 1821, "chrUn_JTFH01000402v1_decoy": 1815, "chrUn_JTFH01000403v1_decoy": 1811, "chrUn_JTFH01000404v1_decoy": 1808, "chrUn_JTFH01000405v1_decoy": 1808, "chrUn_JTFH01000406v1_decoy": 1807, "chrUn_JTFH01000407v1_decoy": 1807, "chrUn_JTFH01000408v1_decoy": 1802, "chrUn_JTFH01000409v1_decoy": 1801, "chrUn_JTFH01000410v1_decoy": 1800, "chrUn_JTFH01000411v1_decoy": 1795, "chrUn_JTFH01000412v1_decoy": 1794, "chrUn_JTFH01000413v1_decoy": 1792, "chrUn_JTFH01000415v1_decoy": 1786, "chrUn_JTFH01000416v1_decoy": 1782, "chrUn_JTFH01000417v1_decoy": 1782, "chrUn_JTFH01000418v1_decoy": 1781, "chrUn_JTFH01000419v1_decoy": 1781, "chrUn_JTFH01000420v1_decoy": 1779, "chrUn_JTFH01000421v1_decoy": 1777, "chrUn_JTFH01000422v1_decoy": 1764, "chrUn_JTFH01000423v1_decoy": 1762, "chrUn_JTFH01000424v1_decoy": 1755, "chrUn_JTFH01000425v1_decoy": 1749, "chrUn_JTFH01000426v1_decoy": 1747, "chrUn_JTFH01000427v1_decoy": 1746, "chrUn_JTFH01000428v1_decoy": 1745, "chrUn_JTFH01000429v1_decoy": 1744, "chrUn_JTFH01000430v1_decoy": 1742, "chrUn_JTFH01000431v1_decoy": 1740, "chrUn_JTFH01000432v1_decoy": 1740, "chrUn_JTFH01000433v1_decoy": 1736, "chrUn_JTFH01000434v1_decoy": 1735, "chrUn_JTFH01000435v1_decoy": 1732, "chrUn_JTFH01000436v1_decoy": 1732, "chrUn_JTFH01000437v1_decoy": 1730, "chrUn_JTFH01000438v1_decoy": 1727, "chrUn_JTFH01000439v1_decoy": 1722, "chrUn_JTFH01000440v1_decoy": 1718, "chrUn_JTFH01000441v1_decoy": 1716, "chrUn_JTFH01000442v1_decoy": 1710, "chrUn_JTFH01000443v1_decoy": 1708, "chrUn_JTFH01000444v1_decoy": 1707, "chrUn_JTFH01000445v1_decoy": 1706, "chrUn_JTFH01000446v1_decoy": 1705, "chrUn_JTFH01000447v1_decoy": 1704, "chrUn_JTFH01000448v1_decoy": 1699, "chrUn_JTFH01000449v1_decoy": 1698, "chrUn_JTFH01000450v1_decoy": 1697, "chrUn_JTFH01000451v1_decoy": 1697, "chrUn_JTFH01000452v1_decoy": 1695, "chrUn_JTFH01000453v1_decoy": 1695, "chrUn_JTFH01000454v1_decoy": 1693, "chrUn_JTFH01000455v1_decoy": 1687, "chrUn_JTFH01000456v1_decoy": 1686, "chrUn_JTFH01000457v1_decoy": 1680, "chrUn_JTFH01000458v1_decoy": 1679, "chrUn_JTFH01000459v1_decoy": 1679, "chrUn_JTFH01000460v1_decoy": 1678, "chrUn_JTFH01000461v1_decoy": 1674, "chrUn_JTFH01000462v1_decoy": 1674, "chrUn_JTFH01000463v1_decoy": 1671, "chrUn_JTFH01000464v1_decoy": 1669, "chrUn_JTFH01000465v1_decoy": 1665, "chrUn_JTFH01000466v1_decoy": 1663, "chrUn_JTFH01000467v1_decoy": 1657, "chrUn_JTFH01000468v1_decoy": 1653, "chrUn_JTFH01000471v1_decoy": 1649, "chrUn_JTFH01000472v1_decoy": 1649, "chrUn_JTFH01000473v1_decoy": 1640, "chrUn_JTFH01000474v1_decoy": 1638, "chrUn_JTFH01000475v1_decoy": 1636, "chrUn_JTFH01000476v1_decoy": 1632, "chrUn_JTFH01000477v1_decoy": 1631, "chrUn_JTFH01000478v1_decoy": 1630, "chrUn_JTFH01000479v1_decoy": 1627, "chrUn_JTFH01000480v1_decoy": 1624, "chrUn_JTFH01000481v1_decoy": 1617, "chrUn_JTFH01000482v1_decoy": 1616, "chrUn_JTFH01000483v1_decoy": 1615, "chrUn_JTFH01000484v1_decoy": 1611, "chrUn_JTFH01000485v1_decoy": 1611, "chrUn_JTFH01000486v1_decoy": 1606, "chrUn_JTFH01000487v1_decoy": 1605, "chrUn_JTFH01000488v1_decoy": 1605, "chrUn_JTFH01000489v1_decoy": 1600, "chrUn_JTFH01000490v1_decoy": 1598, "chrUn_JTFH01000491v1_decoy": 1598, "chrUn_JTFH01000492v1_decoy": 1597, "chrUn_JTFH01000493v1_decoy": 1596, "chrUn_JTFH01000494v1_decoy": 1595, "chrUn_JTFH01000495v1_decoy": 1592, "chrUn_JTFH01000496v1_decoy": 1589, "chrUn_JTFH01000497v1_decoy": 1585, "chrUn_JTFH01000498v1_decoy": 1579, "chrUn_JTFH01000499v1_decoy": 1578, "chrUn_JTFH01000500v1_decoy": 1577, "chrUn_JTFH01000501v1_decoy": 1577, "chrUn_JTFH01000502v1_decoy": 1577, "chrUn_JTFH01000503v1_decoy": 1576, "chrUn_JTFH01000504v1_decoy": 1575, "chrUn_JTFH01000505v1_decoy": 1574, "chrUn_JTFH01000506v1_decoy": 1572, "chrUn_JTFH01000507v1_decoy": 1571, "chrUn_JTFH01000508v1_decoy": 1563, "chrUn_JTFH01000509v1_decoy": 1561, "chrUn_JTFH01000510v1_decoy": 1561, "chrUn_JTFH01000511v1_decoy": 1560, "chrUn_JTFH01000512v1_decoy": 1560, "chrUn_JTFH01000513v1_decoy": 1554, "chrUn_JTFH01000514v1_decoy": 1552, "chrUn_JTFH01000515v1_decoy": 1548, "chrUn_JTFH01000516v1_decoy": 1546, "chrUn_JTFH01000517v1_decoy": 1541, "chrUn_JTFH01000518v1_decoy": 1536, "chrUn_JTFH01000519v1_decoy": 1533, "chrUn_JTFH01000520v1_decoy": 1532, "chrUn_JTFH01000521v1_decoy": 1532, "chrUn_JTFH01000522v1_decoy": 1530, "chrUn_JTFH01000523v1_decoy": 1527, "chrUn_JTFH01000524v1_decoy": 1526, "chrUn_JTFH01000525v1_decoy": 1524, "chrUn_JTFH01000526v1_decoy": 1523, "chrUn_JTFH01000527v1_decoy": 1523, "chrUn_JTFH01000528v1_decoy": 1522, "chrUn_JTFH01000529v1_decoy": 1522, "chrUn_JTFH01000530v1_decoy": 1519, "chrUn_JTFH01000531v1_decoy": 1513, "chrUn_JTFH01000532v1_decoy": 1508, "chrUn_JTFH01000533v1_decoy": 1508, "chrUn_JTFH01000534v1_decoy": 1505, "chrUn_JTFH01000535v1_decoy": 1503, "chrUn_JTFH01000536v1_decoy": 1496, "chrUn_JTFH01000537v1_decoy": 1491, "chrUn_JTFH01000538v1_decoy": 1490, "chrUn_JTFH01000539v1_decoy": 1490, "chrUn_JTFH01000540v1_decoy": 1487, "chrUn_JTFH01000541v1_decoy": 1486, "chrUn_JTFH01000542v1_decoy": 1485, "chrUn_JTFH01000544v1_decoy": 1483, "chrUn_JTFH01000545v1_decoy": 1479, "chrUn_JTFH01000546v1_decoy": 1479, "chrUn_JTFH01000547v1_decoy": 1476, "chrUn_JTFH01000548v1_decoy": 1475, "chrUn_JTFH01000550v1_decoy": 1469, "chrUn_JTFH01000551v1_decoy": 1468, "chrUn_JTFH01000552v1_decoy": 1467, "chrUn_JTFH01000553v1_decoy": 1465, "chrUn_JTFH01000554v1_decoy": 1464, "chrUn_JTFH01000555v1_decoy": 1463, "chrUn_JTFH01000556v1_decoy": 1463, "chrUn_JTFH01000557v1_decoy": 1459, "chrUn_JTFH01000558v1_decoy": 1459, "chrUn_JTFH01000559v1_decoy": 1458, "chrUn_JTFH01000560v1_decoy": 1458, "chrUn_JTFH01000561v1_decoy": 1454, "chrUn_JTFH01000562v1_decoy": 1449, "chrUn_JTFH01000563v1_decoy": 1449, "chrUn_JTFH01000564v1_decoy": 1448, "chrUn_JTFH01000565v1_decoy": 1446, "chrUn_JTFH01000566v1_decoy": 1442, "chrUn_JTFH01000567v1_decoy": 1441, "chrUn_JTFH01000568v1_decoy": 1440, "chrUn_JTFH01000569v1_decoy": 1439, "chrUn_JTFH01000570v1_decoy": 1437, "chrUn_JTFH01000571v1_decoy": 1436, "chrUn_JTFH01000572v1_decoy": 1429, "chrUn_JTFH01000573v1_decoy": 1429, "chrUn_JTFH01000574v1_decoy": 1427, "chrUn_JTFH01000575v1_decoy": 1426, "chrUn_JTFH01000576v1_decoy": 1425, "chrUn_JTFH01000577v1_decoy": 1424, "chrUn_JTFH01000578v1_decoy": 1424, "chrUn_JTFH01000579v1_decoy": 1423, "chrUn_JTFH01000580v1_decoy": 1423, "chrUn_JTFH01000581v1_decoy": 1423, "chrUn_JTFH01000582v1_decoy": 1414, "chrUn_JTFH01000583v1_decoy": 1414, "chrUn_JTFH01000584v1_decoy": 1413, "chrUn_JTFH01000585v1_decoy": 1413, "chrUn_JTFH01000586v1_decoy": 1410, "chrUn_JTFH01000587v1_decoy": 1409, "chrUn_JTFH01000588v1_decoy": 1409, "chrUn_JTFH01000589v1_decoy": 1406, "chrUn_JTFH01000590v1_decoy": 1405, "chrUn_JTFH01000591v1_decoy": 1405, "chrUn_JTFH01000592v1_decoy": 1404, "chrUn_JTFH01000593v1_decoy": 1404, "chrUn_JTFH01000594v1_decoy": 1402, "chrUn_JTFH01000595v1_decoy": 1402, "chrUn_JTFH01000596v1_decoy": 1402, "chrUn_JTFH01000597v1_decoy": 1402, "chrUn_JTFH01000599v1_decoy": 1398, "chrUn_JTFH01000600v1_decoy": 1396, "chrUn_JTFH01000601v1_decoy": 1395, "chrUn_JTFH01000602v1_decoy": 1394, "chrUn_JTFH01000603v1_decoy": 1393, "chrUn_JTFH01000604v1_decoy": 1391, "chrUn_JTFH01000605v1_decoy": 1389, "chrUn_JTFH01000606v1_decoy": 1389, "chrUn_JTFH01000607v1_decoy": 1388, "chrUn_JTFH01000608v1_decoy": 1387, "chrUn_JTFH01000609v1_decoy": 1384, "chrUn_JTFH01000610v1_decoy": 1381, "chrUn_JTFH01000611v1_decoy": 1381, "chrUn_JTFH01000612v1_decoy": 1379, "chrUn_JTFH01000613v1_decoy": 1377, "chrUn_JTFH01000614v1_decoy": 1376, "chrUn_JTFH01000615v1_decoy": 1376, "chrUn_JTFH01000616v1_decoy": 1375, "chrUn_JTFH01000617v1_decoy": 1374, "chrUn_JTFH01000618v1_decoy": 1372, "chrUn_JTFH01000619v1_decoy": 1371, "chrUn_JTFH01000620v1_decoy": 1370, "chrUn_JTFH01000621v1_decoy": 1370, "chrUn_JTFH01000622v1_decoy": 1366, "chrUn_JTFH01000623v1_decoy": 1363, "chrUn_JTFH01000624v1_decoy": 1360, "chrUn_JTFH01000625v1_decoy": 1356, "chrUn_JTFH01000626v1_decoy": 1355, "chrUn_JTFH01000627v1_decoy": 1355, "chrUn_JTFH01000628v1_decoy": 1352, "chrUn_JTFH01000629v1_decoy": 1345, "chrUn_JTFH01000630v1_decoy": 1344, "chrUn_JTFH01000631v1_decoy": 1344, "chrUn_JTFH01000632v1_decoy": 1342, "chrUn_JTFH01000633v1_decoy": 1342, "chrUn_JTFH01000634v1_decoy": 1336, "chrUn_JTFH01000635v1_decoy": 1334, "chrUn_JTFH01000636v1_decoy": 1334, "chrUn_JTFH01000637v1_decoy": 1333, "chrUn_JTFH01000638v1_decoy": 1332, "chrUn_JTFH01000639v1_decoy": 1328, "chrUn_JTFH01000640v1_decoy": 1328, "chrUn_JTFH01000641v1_decoy": 1328, "chrUn_JTFH01000642v1_decoy": 1327, "chrUn_JTFH01000643v1_decoy": 1325, "chrUn_JTFH01000644v1_decoy": 1322, "chrUn_JTFH01000645v1_decoy": 1320, "chrUn_JTFH01000646v1_decoy": 1319, "chrUn_JTFH01000647v1_decoy": 1318, "chrUn_JTFH01000648v1_decoy": 1315, "chrUn_JTFH01000649v1_decoy": 1314, "chrUn_JTFH01000650v1_decoy": 1313, "chrUn_JTFH01000651v1_decoy": 1313, "chrUn_JTFH01000652v1_decoy": 1312, "chrUn_JTFH01000653v1_decoy": 1310, "chrUn_JTFH01000654v1_decoy": 1309, "chrUn_JTFH01000655v1_decoy": 1309, "chrUn_JTFH01000656v1_decoy": 1307, "chrUn_JTFH01000657v1_decoy": 1307, "chrUn_JTFH01000658v1_decoy": 1305, "chrUn_JTFH01000659v1_decoy": 1304, "chrUn_JTFH01000660v1_decoy": 1303, "chrUn_JTFH01000661v1_decoy": 1302, "chrUn_JTFH01000662v1_decoy": 1302, "chrUn_JTFH01000663v1_decoy": 1301, "chrUn_JTFH01000664v1_decoy": 1301, "chrUn_JTFH01000666v1_decoy": 1299, "chrUn_JTFH01000667v1_decoy": 1297, "chrUn_JTFH01000668v1_decoy": 1295, "chrUn_JTFH01000669v1_decoy": 1294, "chrUn_JTFH01000670v1_decoy": 1293, "chrUn_JTFH01000671v1_decoy": 1291, "chrUn_JTFH01000672v1_decoy": 1291, "chrUn_JTFH01000673v1_decoy": 1289, "chrUn_JTFH01000674v1_decoy": 1288, "chrUn_JTFH01000675v1_decoy": 1288, "chrUn_JTFH01000676v1_decoy": 1287, "chrUn_JTFH01000677v1_decoy": 1287, "chrUn_JTFH01000678v1_decoy": 1287, "chrUn_JTFH01000679v1_decoy": 1286, "chrUn_JTFH01000680v1_decoy": 1283, "chrUn_JTFH01000681v1_decoy": 1281, "chrUn_JTFH01000682v1_decoy": 1277, "chrUn_JTFH01000683v1_decoy": 1274, "chrUn_JTFH01000684v1_decoy": 1270, "chrUn_JTFH01000685v1_decoy": 1267, "chrUn_JTFH01000686v1_decoy": 1266, "chrUn_JTFH01000687v1_decoy": 1260, "chrUn_JTFH01000688v1_decoy": 1259, "chrUn_JTFH01000689v1_decoy": 1258, "chrUn_JTFH01000690v1_decoy": 1258, "chrUn_JTFH01000691v1_decoy": 1258, "chrUn_JTFH01000692v1_decoy": 1256, "chrUn_JTFH01000693v1_decoy": 1255, "chrUn_JTFH01000694v1_decoy": 1254, "chrUn_JTFH01000695v1_decoy": 1254, "chrUn_JTFH01000696v1_decoy": 1253, "chrUn_JTFH01000697v1_decoy": 1250, "chrUn_JTFH01000698v1_decoy": 1249, "chrUn_JTFH01000699v1_decoy": 1248, "chrUn_JTFH01000700v1_decoy": 1248, "chrUn_JTFH01000701v1_decoy": 1247, "chrUn_JTFH01000702v1_decoy": 1242, "chrUn_JTFH01000703v1_decoy": 1242, "chrUn_JTFH01000704v1_decoy": 1241, "chrUn_JTFH01000705v1_decoy": 1241, "chrUn_JTFH01000706v1_decoy": 1241, "chrUn_JTFH01000707v1_decoy": 1239, "chrUn_JTFH01000708v1_decoy": 1238, "chrUn_JTFH01000709v1_decoy": 1237, "chrUn_JTFH01000710v1_decoy": 1236, "chrUn_JTFH01000711v1_decoy": 1235, "chrUn_JTFH01000712v1_decoy": 1234, "chrUn_JTFH01000713v1_decoy": 1234, "chrUn_JTFH01000714v1_decoy": 1234, "chrUn_JTFH01000716v1_decoy": 1232, "chrUn_JTFH01000717v1_decoy": 1232, "chrUn_JTFH01000718v1_decoy": 1231, "chrUn_JTFH01000719v1_decoy": 1230, "chrUn_JTFH01000720v1_decoy": 1228, "chrUn_JTFH01000721v1_decoy": 1227, "chrUn_JTFH01000722v1_decoy": 1227, "chrUn_JTFH01000723v1_decoy": 1226, "chrUn_JTFH01000724v1_decoy": 1224, "chrUn_JTFH01000725v1_decoy": 1224, "chrUn_JTFH01000726v1_decoy": 1220, "chrUn_JTFH01000727v1_decoy": 1220, "chrUn_JTFH01000728v1_decoy": 1219, "chrUn_JTFH01000729v1_decoy": 1217, "chrUn_JTFH01000731v1_decoy": 1215, "chrUn_JTFH01000732v1_decoy": 1214, "chrUn_JTFH01000733v1_decoy": 1214, "chrUn_JTFH01000734v1_decoy": 1214, "chrUn_JTFH01000735v1_decoy": 1213, "chrUn_JTFH01000736v1_decoy": 1212, "chrUn_JTFH01000737v1_decoy": 1209, "chrUn_JTFH01000738v1_decoy": 1208, "chrUn_JTFH01000739v1_decoy": 1207, "chrUn_JTFH01000740v1_decoy": 1207, "chrUn_JTFH01000741v1_decoy": 1207, "chrUn_JTFH01000742v1_decoy": 1206, "chrUn_JTFH01000743v1_decoy": 1206, "chrUn_JTFH01000744v1_decoy": 1205, "chrUn_JTFH01000745v1_decoy": 1205, "chrUn_JTFH01000746v1_decoy": 1204, "chrUn_JTFH01000747v1_decoy": 1204, "chrUn_JTFH01000748v1_decoy": 1204, "chrUn_JTFH01000749v1_decoy": 1203, "chrUn_JTFH01000752v1_decoy": 1200, "chrUn_JTFH01000753v1_decoy": 1200, "chrUn_JTFH01000754v1_decoy": 1199, "chrUn_JTFH01000755v1_decoy": 1198, "chrUn_JTFH01000756v1_decoy": 1197, "chrUn_JTFH01000757v1_decoy": 1196, "chrUn_JTFH01000758v1_decoy": 1195, "chrUn_JTFH01000759v1_decoy": 1194, "chrUn_JTFH01000760v1_decoy": 1194, "chrUn_JTFH01000761v1_decoy": 1191, "chrUn_JTFH01000762v1_decoy": 1189, "chrUn_JTFH01000763v1_decoy": 1186, "chrUn_JTFH01000764v1_decoy": 1186, "chrUn_JTFH01000765v1_decoy": 1184, "chrUn_JTFH01000766v1_decoy": 1183, "chrUn_JTFH01000767v1_decoy": 1183, "chrUn_JTFH01000768v1_decoy": 1182, "chrUn_JTFH01000769v1_decoy": 1181, "chrUn_JTFH01000770v1_decoy": 1181, "chrUn_JTFH01000771v1_decoy": 1181, "chrUn_JTFH01000772v1_decoy": 1181, "chrUn_JTFH01000774v1_decoy": 1178, "chrUn_JTFH01000775v1_decoy": 1178, "chrUn_JTFH01000776v1_decoy": 1177, "chrUn_JTFH01000777v1_decoy": 1177, "chrUn_JTFH01000778v1_decoy": 1171, "chrUn_JTFH01000779v1_decoy": 1171, "chrUn_JTFH01000780v1_decoy": 1171, "chrUn_JTFH01000781v1_decoy": 1170, "chrUn_JTFH01000782v1_decoy": 1170, "chrUn_JTFH01000783v1_decoy": 1167, "chrUn_JTFH01000784v1_decoy": 1167, "chrUn_JTFH01000785v1_decoy": 1167, "chrUn_JTFH01000786v1_decoy": 1165, "chrUn_JTFH01000787v1_decoy": 1165, "chrUn_JTFH01000788v1_decoy": 1162, "chrUn_JTFH01000789v1_decoy": 1157, "chrUn_JTFH01000790v1_decoy": 1156, "chrUn_JTFH01000791v1_decoy": 1156, "chrUn_JTFH01000792v1_decoy": 1154, "chrUn_JTFH01000793v1_decoy": 1154, "chrUn_JTFH01000794v1_decoy": 1151, "chrUn_JTFH01000795v1_decoy": 1151, "chrUn_JTFH01000796v1_decoy": 1150, "chrUn_JTFH01000797v1_decoy": 1150, "chrUn_JTFH01000798v1_decoy": 1147, "chrUn_JTFH01000799v1_decoy": 1147, "chrUn_JTFH01000800v1_decoy": 1146, "chrUn_JTFH01000801v1_decoy": 1144, "chrUn_JTFH01000802v1_decoy": 1144, "chrUn_JTFH01000804v1_decoy": 1142, "chrUn_JTFH01000805v1_decoy": 1141, "chrUn_JTFH01000806v1_decoy": 1141, "chrUn_JTFH01000807v1_decoy": 1140, "chrUn_JTFH01000808v1_decoy": 1138, "chrUn_JTFH01000809v1_decoy": 1134, "chrUn_JTFH01000810v1_decoy": 1134, "chrUn_JTFH01000811v1_decoy": 1132, "chrUn_JTFH01000812v1_decoy": 1131, "chrUn_JTFH01000813v1_decoy": 1131, "chrUn_JTFH01000814v1_decoy": 1130, "chrUn_JTFH01000815v1_decoy": 1127, "chrUn_JTFH01000816v1_decoy": 1126, "chrUn_JTFH01000817v1_decoy": 1124, "chrUn_JTFH01000818v1_decoy": 1122, "chrUn_JTFH01000819v1_decoy": 1122, "chrUn_JTFH01000821v1_decoy": 1119, "chrUn_JTFH01000822v1_decoy": 1119, "chrUn_JTFH01000823v1_decoy": 1119, "chrUn_JTFH01000824v1_decoy": 1119, "chrUn_JTFH01000825v1_decoy": 1118, "chrUn_JTFH01000826v1_decoy": 1116, "chrUn_JTFH01000827v1_decoy": 1116, "chrUn_JTFH01000828v1_decoy": 1115, "chrUn_JTFH01000829v1_decoy": 1115, "chrUn_JTFH01000830v1_decoy": 1115, "chrUn_JTFH01000831v1_decoy": 1114, "chrUn_JTFH01000832v1_decoy": 1113, "chrUn_JTFH01000833v1_decoy": 1113, "chrUn_JTFH01000834v1_decoy": 1110, "chrUn_JTFH01000835v1_decoy": 1110, "chrUn_JTFH01000836v1_decoy": 1109, "chrUn_JTFH01000837v1_decoy": 1108, "chrUn_JTFH01000838v1_decoy": 1107, "chrUn_JTFH01000839v1_decoy": 1107, "chrUn_JTFH01000840v1_decoy": 1107, "chrUn_JTFH01000841v1_decoy": 1107, "chrUn_JTFH01000842v1_decoy": 1106, "chrUn_JTFH01000843v1_decoy": 1103, "chrUn_JTFH01000844v1_decoy": 1103, "chrUn_JTFH01000845v1_decoy": 1103, "chrUn_JTFH01000846v1_decoy": 1100, "chrUn_JTFH01000847v1_decoy": 1099, "chrUn_JTFH01000848v1_decoy": 1098, "chrUn_JTFH01000849v1_decoy": 1097, "chrUn_JTFH01000850v1_decoy": 1096, "chrUn_JTFH01000851v1_decoy": 1096, "chrUn_JTFH01000852v1_decoy": 1094, "chrUn_JTFH01000853v1_decoy": 1093, "chrUn_JTFH01000854v1_decoy": 1090, "chrUn_JTFH01000855v1_decoy": 1088, "chrUn_JTFH01000856v1_decoy": 1087, "chrUn_JTFH01000857v1_decoy": 1086, "chrUn_JTFH01000858v1_decoy": 1085, "chrUn_JTFH01000859v1_decoy": 1084, "chrUn_JTFH01000860v1_decoy": 1084, "chrUn_JTFH01000861v1_decoy": 1084, "chrUn_JTFH01000862v1_decoy": 1084, "chrUn_JTFH01000863v1_decoy": 1083, "chrUn_JTFH01000864v1_decoy": 1083, "chrUn_JTFH01000865v1_decoy": 1082, "chrUn_JTFH01000866v1_decoy": 1082, "chrUn_JTFH01000867v1_decoy": 1081, "chrUn_JTFH01000868v1_decoy": 1081, "chrUn_JTFH01000869v1_decoy": 1079, "chrUn_JTFH01000870v1_decoy": 1076, "chrUn_JTFH01000871v1_decoy": 1074, "chrUn_JTFH01000872v1_decoy": 1073, "chrUn_JTFH01000873v1_decoy": 1073, "chrUn_JTFH01000874v1_decoy": 1071, "chrUn_JTFH01000875v1_decoy": 1069, "chrUn_JTFH01000876v1_decoy": 1067, "chrUn_JTFH01000877v1_decoy": 1067, "chrUn_JTFH01000878v1_decoy": 1067, "chrUn_JTFH01000879v1_decoy": 1066, "chrUn_JTFH01000880v1_decoy": 1065, "chrUn_JTFH01000881v1_decoy": 1065, "chrUn_JTFH01000882v1_decoy": 1065, "chrUn_JTFH01000883v1_decoy": 1065, "chrUn_JTFH01000884v1_decoy": 1065, "chrUn_JTFH01000885v1_decoy": 1064, "chrUn_JTFH01000886v1_decoy": 1064, "chrUn_JTFH01000887v1_decoy": 1064, "chrUn_JTFH01000888v1_decoy": 1063, "chrUn_JTFH01000889v1_decoy": 1062, "chrUn_JTFH01000890v1_decoy": 1062, "chrUn_JTFH01000891v1_decoy": 1062, "chrUn_JTFH01000892v1_decoy": 1061, "chrUn_JTFH01000893v1_decoy": 1060, "chrUn_JTFH01000894v1_decoy": 1057, "chrUn_JTFH01000895v1_decoy": 1057, "chrUn_JTFH01000896v1_decoy": 1056, "chrUn_JTFH01000897v1_decoy": 1055, "chrUn_JTFH01000898v1_decoy": 1055, "chrUn_JTFH01000899v1_decoy": 1055, "chrUn_JTFH01000900v1_decoy": 1055, "chrUn_JTFH01000901v1_decoy": 1054, "chrUn_JTFH01000902v1_decoy": 1051, "chrUn_JTFH01000903v1_decoy": 1050, "chrUn_JTFH01000904v1_decoy": 1050, "chrUn_JTFH01000905v1_decoy": 1049, "chrUn_JTFH01000907v1_decoy": 1047, "chrUn_JTFH01000908v1_decoy": 1046, "chrUn_JTFH01000909v1_decoy": 1046, "chrUn_JTFH01000910v1_decoy": 1046, "chrUn_JTFH01000914v1_decoy": 1044, "chrUn_JTFH01000915v1_decoy": 1042, "chrUn_JTFH01000916v1_decoy": 1041, "chrUn_JTFH01000917v1_decoy": 1039, "chrUn_JTFH01000918v1_decoy": 1039, "chrUn_JTFH01000919v1_decoy": 1038, "chrUn_JTFH01000920v1_decoy": 1036, "chrUn_JTFH01000921v1_decoy": 1036, "chrUn_JTFH01000922v1_decoy": 1035, "chrUn_JTFH01000923v1_decoy": 1035, "chrUn_JTFH01000924v1_decoy": 1033, "chrUn_JTFH01000925v1_decoy": 1032, "chrUn_JTFH01000926v1_decoy": 1031, "chrUn_JTFH01000927v1_decoy": 1031, "chrUn_JTFH01000928v1_decoy": 1031, "chrUn_JTFH01000929v1_decoy": 1027, "chrUn_JTFH01000930v1_decoy": 1027, "chrUn_JTFH01000933v1_decoy": 1024, "chrUn_JTFH01000934v1_decoy": 1024, "chrUn_JTFH01000935v1_decoy": 1022, "chrUn_JTFH01000936v1_decoy": 1022, "chrUn_JTFH01000937v1_decoy": 1021, "chrUn_JTFH01000938v1_decoy": 1020, "chrUn_JTFH01000939v1_decoy": 1019, "chrUn_JTFH01000940v1_decoy": 1018, "chrUn_JTFH01000941v1_decoy": 1018, "chrUn_JTFH01000942v1_decoy": 1018, "chrUn_JTFH01000943v1_decoy": 1016, "chrUn_JTFH01000944v1_decoy": 1010, "chrUn_JTFH01000945v1_decoy": 1010, "chrUn_JTFH01000946v1_decoy": 1009, "chrUn_JTFH01000947v1_decoy": 1008, "chrUn_JTFH01000948v1_decoy": 1007, "chrUn_JTFH01000949v1_decoy": 1006, "chrUn_JTFH01000950v1_decoy": 1005, "chrUn_JTFH01000951v1_decoy": 1005, "chrUn_JTFH01000952v1_decoy": 1004, "chrUn_JTFH01000953v1_decoy": 1004, "chrUn_JTFH01000954v1_decoy": 1003, "chrUn_JTFH01000955v1_decoy": 1003, "chrUn_JTFH01000956v1_decoy": 1003, "chrUn_JTFH01000957v1_decoy": 1003, "chrUn_JTFH01000958v1_decoy": 1002, "chrUn_JTFH01000959v1_decoy": 1002, "chrUn_JTFH01000960v1_decoy": 1000, "chrUn_JTFH01000961v1_decoy": 1000, "chrUn_JTFH01000962v1_decoy": 8358, "chrUn_JTFH01000963v1_decoy": 7932, "chrUn_JTFH01000964v1_decoy": 6846, "chrUn_JTFH01000965v1_decoy": 4591, "chrUn_JTFH01000966v1_decoy": 4041, "chrUn_JTFH01000967v1_decoy": 3841, "chrUn_JTFH01000968v1_decoy": 3754, "chrUn_JTFH01000969v1_decoy": 3743, "chrUn_JTFH01000970v1_decoy": 3702, "chrUn_JTFH01000971v1_decoy": 3625, "chrUn_JTFH01000972v1_decoy": 3529, "chrUn_JTFH01000973v1_decoy": 3508, "chrUn_JTFH01000974v1_decoy": 3359, "chrUn_JTFH01000975v1_decoy": 3320, "chrUn_JTFH01000976v1_decoy": 3231, "chrUn_JTFH01000977v1_decoy": 3220, "chrUn_JTFH01000978v1_decoy": 3212, "chrUn_JTFH01000979v1_decoy": 3192, "chrUn_JTFH01000980v1_decoy": 3092, "chrUn_JTFH01000981v1_decoy": 3087, "chrUn_JTFH01000982v1_decoy": 3048, "chrUn_JTFH01000983v1_decoy": 3005, "chrUn_JTFH01000984v1_decoy": 3004, "chrUn_JTFH01000985v1_decoy": 2959, "chrUn_JTFH01000986v1_decoy": 2934, "chrUn_JTFH01000987v1_decoy": 2933, "chrUn_JTFH01000988v1_decoy": 2827, "chrUn_JTFH01000989v1_decoy": 2794, "chrUn_JTFH01000990v1_decoy": 2749, "chrUn_JTFH01000991v1_decoy": 2745, "chrUn_JTFH01000992v1_decoy": 2733, "chrUn_JTFH01000993v1_decoy": 2698, "chrUn_JTFH01000994v1_decoy": 2665, "chrUn_JTFH01000995v1_decoy": 2634, "chrUn_JTFH01000996v1_decoy": 2492, "chrUn_JTFH01000998v1_decoy": 2468, "chrUn_JTFH01000999v1_decoy": 2414, "chrUn_JTFH01001000v1_decoy": 2395, "chrUn_JTFH01001001v1_decoy": 2356, "chrUn_JTFH01001002v1_decoy": 2339, "chrUn_JTFH01001003v1_decoy": 2310, "chrUn_JTFH01001004v1_decoy": 2288, "chrUn_JTFH01001005v1_decoy": 2285, "chrUn_JTFH01001006v1_decoy": 2269, "chrUn_JTFH01001007v1_decoy": 2253, "chrUn_JTFH01001008v1_decoy": 2203, "chrUn_JTFH01001009v1_decoy": 2176, "chrUn_JTFH01001010v1_decoy": 2159, "chrUn_JTFH01001011v1_decoy": 2155, "chrUn_JTFH01001012v1_decoy": 2149, "chrUn_JTFH01001013v1_decoy": 2129, "chrUn_JTFH01001014v1_decoy": 2116, "chrUn_JTFH01001015v1_decoy": 2113, "chrUn_JTFH01001016v1_decoy": 2098, "chrUn_JTFH01001017v1_decoy": 2066, "chrUn_JTFH01001018v1_decoy": 2066, "chrUn_JTFH01001019v1_decoy": 2059, "chrUn_JTFH01001020v1_decoy": 2047, "chrUn_JTFH01001021v1_decoy": 2040, "chrUn_JTFH01001022v1_decoy": 2030, "chrUn_JTFH01001023v1_decoy": 2024, "chrUn_JTFH01001024v1_decoy": 2001, "chrUn_JTFH01001025v1_decoy": 1992, "chrUn_JTFH01001026v1_decoy": 1981, "chrUn_JTFH01001027v1_decoy": 1979, "chrUn_JTFH01001028v1_decoy": 1957, "chrUn_JTFH01001029v1_decoy": 1953, "chrUn_JTFH01001030v1_decoy": 1944, "chrUn_JTFH01001031v1_decoy": 1936, "chrUn_JTFH01001032v1_decoy": 1932, "chrUn_JTFH01001033v1_decoy": 1882, "chrUn_JTFH01001034v1_decoy": 1878, "chrUn_JTFH01001035v1_decoy": 1870, "chrUn_JTFH01001036v1_decoy": 1821, "chrUn_JTFH01001037v1_decoy": 1813, "chrUn_JTFH01001038v1_decoy": 1809, "chrUn_JTFH01001039v1_decoy": 1804, "chrUn_JTFH01001040v1_decoy": 1797, "chrUn_JTFH01001041v1_decoy": 1791, "chrUn_JTFH01001042v1_decoy": 1781, "chrUn_JTFH01001043v1_decoy": 1766, "chrUn_JTFH01001044v1_decoy": 1764, "chrUn_JTFH01001045v1_decoy": 1743, "chrUn_JTFH01001046v1_decoy": 1741, "chrUn_JTFH01001047v1_decoy": 1709, "chrUn_JTFH01001048v1_decoy": 1706, "chrUn_JTFH01001049v1_decoy": 1701, "chrUn_JTFH01001050v1_decoy": 1689, "chrUn_JTFH01001051v1_decoy": 1646, "chrUn_JTFH01001052v1_decoy": 1641, "chrUn_JTFH01001053v1_decoy": 1639, "chrUn_JTFH01001054v1_decoy": 1636, "chrUn_JTFH01001055v1_decoy": 1632, "chrUn_JTFH01001056v1_decoy": 1629, "chrUn_JTFH01001057v1_decoy": 1623, "chrUn_JTFH01001058v1_decoy": 1622, "chrUn_JTFH01001059v1_decoy": 1622, "chrUn_JTFH01001060v1_decoy": 1619, "chrUn_JTFH01001061v1_decoy": 1606, "chrUn_JTFH01001062v1_decoy": 1593, "chrUn_JTFH01001063v1_decoy": 1592, "chrUn_JTFH01001064v1_decoy": 1558, "chrUn_JTFH01001065v1_decoy": 1545, "chrUn_JTFH01001066v1_decoy": 1542, "chrUn_JTFH01001067v1_decoy": 1540, "chrUn_JTFH01001068v1_decoy": 1529, "chrUn_JTFH01001069v1_decoy": 1518, "chrUn_JTFH01001070v1_decoy": 1515, "chrUn_JTFH01001071v1_decoy": 1513, "chrUn_JTFH01001072v1_decoy": 1507, "chrUn_JTFH01001073v1_decoy": 1504, "chrUn_JTFH01001074v1_decoy": 1499, "chrUn_JTFH01001075v1_decoy": 1495, "chrUn_JTFH01001076v1_decoy": 1495, "chrUn_JTFH01001077v1_decoy": 1492, "chrUn_JTFH01001078v1_decoy": 1492, "chrUn_JTFH01001079v1_decoy": 1489, "chrUn_JTFH01001080v1_decoy": 1485, "chrUn_JTFH01001081v1_decoy": 1483, "chrUn_JTFH01001082v1_decoy": 1473, "chrUn_JTFH01001083v1_decoy": 1470, "chrUn_JTFH01001084v1_decoy": 1463, "chrUn_JTFH01001085v1_decoy": 1460, "chrUn_JTFH01001086v1_decoy": 1458, "chrUn_JTFH01001087v1_decoy": 1456, "chrUn_JTFH01001088v1_decoy": 1453, "chrUn_JTFH01001089v1_decoy": 1443, "chrUn_JTFH01001090v1_decoy": 1441, "chrUn_JTFH01001091v1_decoy": 1426, "chrUn_JTFH01001092v1_decoy": 1425, "chrUn_JTFH01001093v1_decoy": 1418, "chrUn_JTFH01001094v1_decoy": 1413, "chrUn_JTFH01001095v1_decoy": 1413, "chrUn_JTFH01001096v1_decoy": 1412, "chrUn_JTFH01001097v1_decoy": 1407, "chrUn_JTFH01001098v1_decoy": 1406, "chrUn_JTFH01001099v1_decoy": 1396, "chrUn_JTFH01001100v1_decoy": 1390, "chrUn_JTFH01001101v1_decoy": 1382, "chrUn_JTFH01001102v1_decoy": 1376, "chrUn_JTFH01001103v1_decoy": 1375, "chrUn_JTFH01001104v1_decoy": 1371, "chrUn_JTFH01001105v1_decoy": 1367, "chrUn_JTFH01001106v1_decoy": 1364, "chrUn_JTFH01001107v1_decoy": 1356, "chrUn_JTFH01001108v1_decoy": 1355, "chrUn_JTFH01001109v1_decoy": 1352, "chrUn_JTFH01001110v1_decoy": 1350, "chrUn_JTFH01001111v1_decoy": 1346, "chrUn_JTFH01001112v1_decoy": 1345, "chrUn_JTFH01001113v1_decoy": 1340, "chrUn_JTFH01001114v1_decoy": 1330, "chrUn_JTFH01001115v1_decoy": 1329, "chrUn_JTFH01001116v1_decoy": 1324, "chrUn_JTFH01001117v1_decoy": 1316, "chrUn_JTFH01001118v1_decoy": 1307, "chrUn_JTFH01001119v1_decoy": 1304, "chrUn_JTFH01001120v1_decoy": 1304, "chrUn_JTFH01001121v1_decoy": 1303, "chrUn_JTFH01001122v1_decoy": 1301, "chrUn_JTFH01001124v1_decoy": 1297, "chrUn_JTFH01001125v1_decoy": 1296, "chrUn_JTFH01001126v1_decoy": 1290, "chrUn_JTFH01001127v1_decoy": 1284, "chrUn_JTFH01001128v1_decoy": 1282, "chrUn_JTFH01001129v1_decoy": 1281, "chrUn_JTFH01001130v1_decoy": 1280, "chrUn_JTFH01001131v1_decoy": 1279, "chrUn_JTFH01001132v1_decoy": 1272, "chrUn_JTFH01001133v1_decoy": 1267, "chrUn_JTFH01001134v1_decoy": 1267, "chrUn_JTFH01001135v1_decoy": 1266, "chrUn_JTFH01001136v1_decoy": 1264, "chrUn_JTFH01001137v1_decoy": 1264, "chrUn_JTFH01001138v1_decoy": 1264, "chrUn_JTFH01001139v1_decoy": 1263, "chrUn_JTFH01001140v1_decoy": 1249, "chrUn_JTFH01001141v1_decoy": 1240, "chrUn_JTFH01001142v1_decoy": 1239, "chrUn_JTFH01001143v1_decoy": 1235, "chrUn_JTFH01001144v1_decoy": 1235, "chrUn_JTFH01001146v1_decoy": 1232, "chrUn_JTFH01001147v1_decoy": 1230, "chrUn_JTFH01001148v1_decoy": 1226, "chrUn_JTFH01001149v1_decoy": 1223, "chrUn_JTFH01001150v1_decoy": 1214, "chrUn_JTFH01001151v1_decoy": 1213, "chrUn_JTFH01001152v1_decoy": 1211, "chrUn_JTFH01001153v1_decoy": 1209, "chrUn_JTFH01001155v1_decoy": 1199, "chrUn_JTFH01001156v1_decoy": 1197, "chrUn_JTFH01001157v1_decoy": 1193, "chrUn_JTFH01001158v1_decoy": 1191, "chrUn_JTFH01001159v1_decoy": 1187, "chrUn_JTFH01001160v1_decoy": 1186, "chrUn_JTFH01001161v1_decoy": 1184, "chrUn_JTFH01001162v1_decoy": 1184, "chrUn_JTFH01001163v1_decoy": 1182, "chrUn_JTFH01001165v1_decoy": 1173, "chrUn_JTFH01001166v1_decoy": 1169, "chrUn_JTFH01001167v1_decoy": 1167, "chrUn_JTFH01001168v1_decoy": 1166, "chrUn_JTFH01001169v1_decoy": 1165, "chrUn_JTFH01001170v1_decoy": 1164, "chrUn_JTFH01001171v1_decoy": 1163, "chrUn_JTFH01001172v1_decoy": 1158, "chrUn_JTFH01001173v1_decoy": 1158, "chrUn_JTFH01001174v1_decoy": 1157, "chrUn_JTFH01001175v1_decoy": 1157, "chrUn_JTFH01001176v1_decoy": 1157, "chrUn_JTFH01001177v1_decoy": 1155, "chrUn_JTFH01001178v1_decoy": 1154, "chrUn_JTFH01001179v1_decoy": 1149, "chrUn_JTFH01001180v1_decoy": 1148, "chrUn_JTFH01001181v1_decoy": 1148, "chrUn_JTFH01001182v1_decoy": 1146, "chrUn_JTFH01001183v1_decoy": 1144, "chrUn_JTFH01001184v1_decoy": 1140, "chrUn_JTFH01001186v1_decoy": 1134, "chrUn_JTFH01001187v1_decoy": 1133, "chrUn_JTFH01001188v1_decoy": 1129, "chrUn_JTFH01001189v1_decoy": 1127, "chrUn_JTFH01001190v1_decoy": 1127, "chrUn_JTFH01001191v1_decoy": 1118, "chrUn_JTFH01001192v1_decoy": 1110, "chrUn_JTFH01001193v1_decoy": 1104, "chrUn_JTFH01001194v1_decoy": 1104, "chrUn_JTFH01001195v1_decoy": 1101, "chrUn_JTFH01001196v1_decoy": 1098, "chrUn_JTFH01001197v1_decoy": 1096, "chrUn_JTFH01001198v1_decoy": 1094, "chrUn_JTFH01001199v1_decoy": 1091, "chrUn_JTFH01001200v1_decoy": 1089, "chrUn_JTFH01001201v1_decoy": 1086, "chrUn_JTFH01001202v1_decoy": 1085, "chrUn_JTFH01001203v1_decoy": 1084, "chrUn_JTFH01001204v1_decoy": 1083, "chrUn_JTFH01001205v1_decoy": 1083, "chrUn_JTFH01001206v1_decoy": 1079, "chrUn_JTFH01001207v1_decoy": 1076, "chrUn_JTFH01001208v1_decoy": 1069, "chrUn_JTFH01001209v1_decoy": 1068, "chrUn_JTFH01001210v1_decoy": 1067, "chrUn_JTFH01001211v1_decoy": 1067, "chrUn_JTFH01001212v1_decoy": 1067, "chrUn_JTFH01001213v1_decoy": 1063, "chrUn_JTFH01001214v1_decoy": 1062, "chrUn_JTFH01001215v1_decoy": 1059, "chrUn_JTFH01001216v1_decoy": 1058, "chrUn_JTFH01001217v1_decoy": 1058, "chrUn_JTFH01001218v1_decoy": 1055, "chrUn_JTFH01001219v1_decoy": 1054, "chrUn_JTFH01001220v1_decoy": 1054, "chrUn_JTFH01001221v1_decoy": 1053, "chrUn_JTFH01001222v1_decoy": 1053, "chrUn_JTFH01001223v1_decoy": 1052, "chrUn_JTFH01001224v1_decoy": 1051, "chrUn_JTFH01001225v1_decoy": 1049, "chrUn_JTFH01001226v1_decoy": 1047, "chrUn_JTFH01001227v1_decoy": 1044, "chrUn_JTFH01001228v1_decoy": 1043, "chrUn_JTFH01001229v1_decoy": 1043, "chrUn_JTFH01001230v1_decoy": 1042, "chrUn_JTFH01001231v1_decoy": 1042, "chrUn_JTFH01001232v1_decoy": 1041, "chrUn_JTFH01001234v1_decoy": 1039, "chrUn_JTFH01001235v1_decoy": 1038, "chrUn_JTFH01001236v1_decoy": 1037, "chrUn_JTFH01001237v1_decoy": 1037, "chrUn_JTFH01001238v1_decoy": 1035, "chrUn_JTFH01001239v1_decoy": 1027, "chrUn_JTFH01001240v1_decoy": 1021, "chrUn_JTFH01001241v1_decoy": 1021, "chrUn_JTFH01001242v1_decoy": 1019, "chrUn_JTFH01001243v1_decoy": 1019, "chrUn_JTFH01001244v1_decoy": 1016, "chrUn_JTFH01001245v1_decoy": 1014, "chrUn_JTFH01001246v1_decoy": 1013, "chrUn_JTFH01001247v1_decoy": 1009, "chrUn_JTFH01001248v1_decoy": 1008, "chrUn_JTFH01001249v1_decoy": 1007, "chrUn_JTFH01001250v1_decoy": 1004, "chrUn_JTFH01001251v1_decoy": 1004, "chrUn_JTFH01001252v1_decoy": 1003, "chrUn_JTFH01001253v1_decoy": 1001, "chrUn_JTFH01001254v1_decoy": 1000, "chrUn_JTFH01001255v1_decoy": 1000, "chrUn_JTFH01001256v1_decoy": 1000, "chrUn_JTFH01001257v1_decoy": 17929, "chrUn_JTFH01001258v1_decoy": 9749, "chrUn_JTFH01001259v1_decoy": 8053, "chrUn_JTFH01001260v1_decoy": 7826, "chrUn_JTFH01001261v1_decoy": 7768, "chrUn_JTFH01001262v1_decoy": 5691, "chrUn_JTFH01001263v1_decoy": 5444, "chrUn_JTFH01001264v1_decoy": 5077, "chrUn_JTFH01001265v1_decoy": 4990, "chrUn_JTFH01001266v1_decoy": 4545, "chrUn_JTFH01001267v1_decoy": 4544, "chrUn_JTFH01001268v1_decoy": 4202, "chrUn_JTFH01001269v1_decoy": 4195, "chrUn_JTFH01001270v1_decoy": 3807, "chrUn_JTFH01001271v1_decoy": 3741, "chrUn_JTFH01001272v1_decoy": 3699, "chrUn_JTFH01001273v1_decoy": 3640, "chrUn_JTFH01001274v1_decoy": 3531, "chrUn_JTFH01001275v1_decoy": 3455, "chrUn_JTFH01001276v1_decoy": 3411, "chrUn_JTFH01001277v1_decoy": 3387, "chrUn_JTFH01001278v1_decoy": 3358, "chrUn_JTFH01001279v1_decoy": 3285, "chrUn_JTFH01001280v1_decoy": 3273, "chrUn_JTFH01001281v1_decoy": 3262, "chrUn_JTFH01001282v1_decoy": 3259, "chrUn_JTFH01001283v1_decoy": 3222, "chrUn_JTFH01001284v1_decoy": 3127, "chrUn_JTFH01001285v1_decoy": 3110, "chrUn_JTFH01001286v1_decoy": 3104, "chrUn_JTFH01001287v1_decoy": 3071, "chrUn_JTFH01001288v1_decoy": 3063, "chrUn_JTFH01001289v1_decoy": 3059, "chrUn_JTFH01001290v1_decoy": 2990, "chrUn_JTFH01001291v1_decoy": 2986, "chrUn_JTFH01001292v1_decoy": 2928, "chrUn_JTFH01001293v1_decoy": 2922, "chrUn_JTFH01001294v1_decoy": 2875, "chrUn_JTFH01001295v1_decoy": 2859, "chrUn_JTFH01001296v1_decoy": 2850, "chrUn_JTFH01001297v1_decoy": 2813, "chrUn_JTFH01001298v1_decoy": 2785, "chrUn_JTFH01001299v1_decoy": 2736, "chrUn_JTFH01001300v1_decoy": 2688, "chrUn_JTFH01001301v1_decoy": 2658, "chrUn_JTFH01001302v1_decoy": 2643, "chrUn_JTFH01001303v1_decoy": 2618, "chrUn_JTFH01001304v1_decoy": 2605, "chrUn_JTFH01001305v1_decoy": 2583, "chrUn_JTFH01001306v1_decoy": 2534, "chrUn_JTFH01001307v1_decoy": 2512, "chrUn_JTFH01001308v1_decoy": 2500, "chrUn_JTFH01001309v1_decoy": 2481, "chrUn_JTFH01001310v1_decoy": 2478, "chrUn_JTFH01001311v1_decoy": 2473, "chrUn_JTFH01001312v1_decoy": 2467, "chrUn_JTFH01001313v1_decoy": 2442, "chrUn_JTFH01001314v1_decoy": 2430, "chrUn_JTFH01001315v1_decoy": 2417, "chrUn_JTFH01001316v1_decoy": 2408, "chrUn_JTFH01001317v1_decoy": 2395, "chrUn_JTFH01001318v1_decoy": 2352, "chrUn_JTFH01001319v1_decoy": 2337, "chrUn_JTFH01001320v1_decoy": 2322, "chrUn_JTFH01001321v1_decoy": 2307, "chrUn_JTFH01001322v1_decoy": 2306, "chrUn_JTFH01001323v1_decoy": 2292, "chrUn_JTFH01001324v1_decoy": 2271, "chrUn_JTFH01001325v1_decoy": 2265, "chrUn_JTFH01001326v1_decoy": 2260, "chrUn_JTFH01001327v1_decoy": 2240, "chrUn_JTFH01001328v1_decoy": 2238, "chrUn_JTFH01001329v1_decoy": 2228, "chrUn_JTFH01001330v1_decoy": 2215, "chrUn_JTFH01001331v1_decoy": 2205, "chrUn_JTFH01001332v1_decoy": 2191, "chrUn_JTFH01001333v1_decoy": 2191, "chrUn_JTFH01001334v1_decoy": 2190, "chrUn_JTFH01001335v1_decoy": 2184, "chrUn_JTFH01001336v1_decoy": 2166, "chrUn_JTFH01001338v1_decoy": 2162, "chrUn_JTFH01001339v1_decoy": 2146, "chrUn_JTFH01001340v1_decoy": 2116, "chrUn_JTFH01001341v1_decoy": 2112, "chrUn_JTFH01001342v1_decoy": 2108, "chrUn_JTFH01001343v1_decoy": 2106, "chrUn_JTFH01001344v1_decoy": 2106, "chrUn_JTFH01001345v1_decoy": 2106, "chrUn_JTFH01001346v1_decoy": 2097, "chrUn_JTFH01001347v1_decoy": 2081, "chrUn_JTFH01001348v1_decoy": 2058, "chrUn_JTFH01001349v1_decoy": 2055, "chrUn_JTFH01001350v1_decoy": 2054, "chrUn_JTFH01001351v1_decoy": 2037, "chrUn_JTFH01001352v1_decoy": 2032, "chrUn_JTFH01001353v1_decoy": 2032, "chrUn_JTFH01001354v1_decoy": 2020, "chrUn_JTFH01001355v1_decoy": 2018, "chrUn_JTFH01001356v1_decoy": 2014, "chrUn_JTFH01001357v1_decoy": 2001, "chrUn_JTFH01001358v1_decoy": 2001, "chrUn_JTFH01001359v1_decoy": 1991, "chrUn_JTFH01001360v1_decoy": 1990, "chrUn_JTFH01001361v1_decoy": 1983, "chrUn_JTFH01001362v1_decoy": 1981, "chrUn_JTFH01001363v1_decoy": 1981, "chrUn_JTFH01001364v1_decoy": 1979, "chrUn_JTFH01001365v1_decoy": 1963, "chrUn_JTFH01001366v1_decoy": 1932, "chrUn_JTFH01001367v1_decoy": 1929, "chrUn_JTFH01001368v1_decoy": 1881, "chrUn_JTFH01001369v1_decoy": 1874, "chrUn_JTFH01001370v1_decoy": 1849, "chrUn_JTFH01001371v1_decoy": 1849, "chrUn_JTFH01001372v1_decoy": 1833, "chrUn_JTFH01001373v1_decoy": 1832, "chrUn_JTFH01001374v1_decoy": 1826, "chrUn_JTFH01001375v1_decoy": 1814, "chrUn_JTFH01001376v1_decoy": 1814, "chrUn_JTFH01001377v1_decoy": 1791, "chrUn_JTFH01001378v1_decoy": 1789, "chrUn_JTFH01001379v1_decoy": 1786, "chrUn_JTFH01001380v1_decoy": 1778, "chrUn_JTFH01001381v1_decoy": 1776, "chrUn_JTFH01001382v1_decoy": 1762, "chrUn_JTFH01001383v1_decoy": 1758, "chrUn_JTFH01001384v1_decoy": 1757, "chrUn_JTFH01001385v1_decoy": 1754, "chrUn_JTFH01001386v1_decoy": 1752, "chrUn_JTFH01001387v1_decoy": 1751, "chrUn_JTFH01001388v1_decoy": 1749, "chrUn_JTFH01001389v1_decoy": 1738, "chrUn_JTFH01001390v1_decoy": 1729, "chrUn_JTFH01001391v1_decoy": 1726, "chrUn_JTFH01001392v1_decoy": 1716, "chrUn_JTFH01001393v1_decoy": 1712, "chrUn_JTFH01001394v1_decoy": 1711, "chrUn_JTFH01001395v1_decoy": 1703, "chrUn_JTFH01001396v1_decoy": 1702, "chrUn_JTFH01001397v1_decoy": 1699, "chrUn_JTFH01001398v1_decoy": 1686, "chrUn_JTFH01001399v1_decoy": 1684, "chrUn_JTFH01001400v1_decoy": 1680, "chrUn_JTFH01001401v1_decoy": 1678, "chrUn_JTFH01001402v1_decoy": 1678, "chrUn_JTFH01001403v1_decoy": 1677, "chrUn_JTFH01001404v1_decoy": 1676, "chrUn_JTFH01001405v1_decoy": 1672, "chrUn_JTFH01001406v1_decoy": 1669, "chrUn_JTFH01001407v1_decoy": 1668, "chrUn_JTFH01001408v1_decoy": 1663, "chrUn_JTFH01001409v1_decoy": 1660, "chrUn_JTFH01001410v1_decoy": 1660, "chrUn_JTFH01001412v1_decoy": 1656, "chrUn_JTFH01001413v1_decoy": 1656, "chrUn_JTFH01001415v1_decoy": 1647, "chrUn_JTFH01001416v1_decoy": 1645, "chrUn_JTFH01001417v1_decoy": 1641, "chrUn_JTFH01001418v1_decoy": 1638, "chrUn_JTFH01001419v1_decoy": 1633, "chrUn_JTFH01001420v1_decoy": 1626, "chrUn_JTFH01001421v1_decoy": 1614, "chrUn_JTFH01001422v1_decoy": 1612, "chrUn_JTFH01001423v1_decoy": 1605, "chrUn_JTFH01001424v1_decoy": 1603, "chrUn_JTFH01001426v1_decoy": 1589, "chrUn_JTFH01001427v1_decoy": 1588, "chrUn_JTFH01001428v1_decoy": 1585, "chrUn_JTFH01001429v1_decoy": 1584, "chrUn_JTFH01001430v1_decoy": 1584, "chrUn_JTFH01001431v1_decoy": 1580, "chrUn_JTFH01001432v1_decoy": 1572, "chrUn_JTFH01001433v1_decoy": 1570, "chrUn_JTFH01001434v1_decoy": 1569, "chrUn_JTFH01001435v1_decoy": 1568, "chrUn_JTFH01001436v1_decoy": 1567, "chrUn_JTFH01001437v1_decoy": 1565, "chrUn_JTFH01001438v1_decoy": 1559, "chrUn_JTFH01001439v1_decoy": 1559, "chrUn_JTFH01001440v1_decoy": 1556, "chrUn_JTFH01001441v1_decoy": 1554, "chrUn_JTFH01001442v1_decoy": 1549, "chrUn_JTFH01001443v1_decoy": 1542, "chrUn_JTFH01001444v1_decoy": 1541, "chrUn_JTFH01001445v1_decoy": 1538, "chrUn_JTFH01001447v1_decoy": 1535, "chrUn_JTFH01001448v1_decoy": 1530, "chrUn_JTFH01001449v1_decoy": 1528, "chrUn_JTFH01001450v1_decoy": 1522, "chrUn_JTFH01001451v1_decoy": 1514, "chrUn_JTFH01001452v1_decoy": 1509, "chrUn_JTFH01001453v1_decoy": 1507, "chrUn_JTFH01001454v1_decoy": 1500, "chrUn_JTFH01001455v1_decoy": 1499, "chrUn_JTFH01001456v1_decoy": 1499, "chrUn_JTFH01001457v1_decoy": 1497, "chrUn_JTFH01001458v1_decoy": 1496, "chrUn_JTFH01001459v1_decoy": 1488, "chrUn_JTFH01001460v1_decoy": 1486, "chrUn_JTFH01001461v1_decoy": 1485, "chrUn_JTFH01001462v1_decoy": 1481, "chrUn_JTFH01001463v1_decoy": 1479, "chrUn_JTFH01001466v1_decoy": 1470, "chrUn_JTFH01001467v1_decoy": 1466, "chrUn_JTFH01001468v1_decoy": 1465, "chrUn_JTFH01001469v1_decoy": 1461, "chrUn_JTFH01001470v1_decoy": 1458, "chrUn_JTFH01001471v1_decoy": 1457, "chrUn_JTFH01001472v1_decoy": 1448, "chrUn_JTFH01001473v1_decoy": 1447, "chrUn_JTFH01001475v1_decoy": 1443, "chrUn_JTFH01001476v1_decoy": 1443, "chrUn_JTFH01001477v1_decoy": 1438, "chrUn_JTFH01001478v1_decoy": 1432, "chrUn_JTFH01001479v1_decoy": 1430, "chrUn_JTFH01001480v1_decoy": 1430, "chrUn_JTFH01001481v1_decoy": 1429, "chrUn_JTFH01001482v1_decoy": 1429, "chrUn_JTFH01001483v1_decoy": 1429, "chrUn_JTFH01001484v1_decoy": 1426, "chrUn_JTFH01001485v1_decoy": 1426, "chrUn_JTFH01001486v1_decoy": 1420, "chrUn_JTFH01001487v1_decoy": 1416, "chrUn_JTFH01001488v1_decoy": 1416, "chrUn_JTFH01001489v1_decoy": 1415, "chrUn_JTFH01001490v1_decoy": 1415, "chrUn_JTFH01001491v1_decoy": 1414, "chrUn_JTFH01001492v1_decoy": 1413, "chrUn_JTFH01001493v1_decoy": 1410, "chrUn_JTFH01001494v1_decoy": 1405, "chrUn_JTFH01001495v1_decoy": 1402, "chrUn_JTFH01001496v1_decoy": 1398, "chrUn_JTFH01001497v1_decoy": 1397, "chrUn_JTFH01001498v1_decoy": 1395, "chrUn_JTFH01001499v1_decoy": 1392, "chrUn_JTFH01001500v1_decoy": 1388, "chrUn_JTFH01001501v1_decoy": 1386, "chrUn_JTFH01001502v1_decoy": 1382, "chrUn_JTFH01001503v1_decoy": 1381, "chrUn_JTFH01001504v1_decoy": 1379, "chrUn_JTFH01001505v1_decoy": 1376, "chrUn_JTFH01001506v1_decoy": 1374, "chrUn_JTFH01001507v1_decoy": 1374, "chrUn_JTFH01001508v1_decoy": 1373, "chrUn_JTFH01001509v1_decoy": 1373, "chrUn_JTFH01001510v1_decoy": 1372, "chrUn_JTFH01001511v1_decoy": 1370, "chrUn_JTFH01001512v1_decoy": 1367, "chrUn_JTFH01001513v1_decoy": 1365, "chrUn_JTFH01001514v1_decoy": 1364, "chrUn_JTFH01001517v1_decoy": 1355, "chrUn_JTFH01001518v1_decoy": 1355, "chrUn_JTFH01001519v1_decoy": 1354, "chrUn_JTFH01001520v1_decoy": 1353, "chrUn_JTFH01001521v1_decoy": 1349, "chrUn_JTFH01001522v1_decoy": 1345, "chrUn_JTFH01001523v1_decoy": 1344, "chrUn_JTFH01001524v1_decoy": 1343, "chrUn_JTFH01001525v1_decoy": 1338, "chrUn_JTFH01001526v1_decoy": 1338, "chrUn_JTFH01001527v1_decoy": 1338, "chrUn_JTFH01001528v1_decoy": 1336, "chrUn_JTFH01001529v1_decoy": 1333, "chrUn_JTFH01001530v1_decoy": 1333, "chrUn_JTFH01001531v1_decoy": 1332, "chrUn_JTFH01001532v1_decoy": 1324, "chrUn_JTFH01001533v1_decoy": 1323, "chrUn_JTFH01001534v1_decoy": 1323, "chrUn_JTFH01001535v1_decoy": 1320, "chrUn_JTFH01001536v1_decoy": 1320, "chrUn_JTFH01001537v1_decoy": 1317, "chrUn_JTFH01001538v1_decoy": 1316, "chrUn_JTFH01001539v1_decoy": 1304, "chrUn_JTFH01001540v1_decoy": 1304, "chrUn_JTFH01001541v1_decoy": 1303, "chrUn_JTFH01001542v1_decoy": 1302, "chrUn_JTFH01001543v1_decoy": 1301, "chrUn_JTFH01001546v1_decoy": 1297, "chrUn_JTFH01001547v1_decoy": 1295, "chrUn_JTFH01001548v1_decoy": 1284, "chrUn_JTFH01001549v1_decoy": 1283, "chrUn_JTFH01001550v1_decoy": 1283, "chrUn_JTFH01001551v1_decoy": 1279, "chrUn_JTFH01001552v1_decoy": 1278, "chrUn_JTFH01001553v1_decoy": 1271, "chrUn_JTFH01001554v1_decoy": 1271, "chrUn_JTFH01001555v1_decoy": 1268, "chrUn_JTFH01001556v1_decoy": 1264, "chrUn_JTFH01001557v1_decoy": 1263, "chrUn_JTFH01001558v1_decoy": 1262, "chrUn_JTFH01001559v1_decoy": 1261, "chrUn_JTFH01001560v1_decoy": 1260, "chrUn_JTFH01001561v1_decoy": 1259, "chrUn_JTFH01001562v1_decoy": 1259, "chrUn_JTFH01001563v1_decoy": 1258, "chrUn_JTFH01001564v1_decoy": 1256, "chrUn_JTFH01001565v1_decoy": 1253, "chrUn_JTFH01001566v1_decoy": 1248, "chrUn_JTFH01001567v1_decoy": 1248, "chrUn_JTFH01001568v1_decoy": 1246, "chrUn_JTFH01001569v1_decoy": 1246, "chrUn_JTFH01001570v1_decoy": 1244, "chrUn_JTFH01001571v1_decoy": 1238, "chrUn_JTFH01001572v1_decoy": 1238, "chrUn_JTFH01001573v1_decoy": 1236, "chrUn_JTFH01001574v1_decoy": 1234, "chrUn_JTFH01001575v1_decoy": 1234, "chrUn_JTFH01001576v1_decoy": 1231, "chrUn_JTFH01001577v1_decoy": 1231, "chrUn_JTFH01001578v1_decoy": 1230, "chrUn_JTFH01001579v1_decoy": 1230, "chrUn_JTFH01001580v1_decoy": 1228, "chrUn_JTFH01001581v1_decoy": 1227, "chrUn_JTFH01001582v1_decoy": 1222, "chrUn_JTFH01001583v1_decoy": 1222, "chrUn_JTFH01001584v1_decoy": 1221, "chrUn_JTFH01001585v1_decoy": 1221, "chrUn_JTFH01001586v1_decoy": 1220, "chrUn_JTFH01001587v1_decoy": 1218, "chrUn_JTFH01001588v1_decoy": 1218, "chrUn_JTFH01001591v1_decoy": 1212, "chrUn_JTFH01001592v1_decoy": 1210, "chrUn_JTFH01001593v1_decoy": 1209, "chrUn_JTFH01001594v1_decoy": 1208, "chrUn_JTFH01001595v1_decoy": 1208, "chrUn_JTFH01001596v1_decoy": 1206, "chrUn_JTFH01001597v1_decoy": 1205, "chrUn_JTFH01001598v1_decoy": 1205, "chrUn_JTFH01001600v1_decoy": 1200, "chrUn_JTFH01001601v1_decoy": 1199, "chrUn_JTFH01001602v1_decoy": 1198, "chrUn_JTFH01001603v1_decoy": 1198, "chrUn_JTFH01001604v1_decoy": 1198, "chrUn_JTFH01001605v1_decoy": 1195, "chrUn_JTFH01001606v1_decoy": 1194, "chrUn_JTFH01001607v1_decoy": 1191, "chrUn_JTFH01001608v1_decoy": 1189, "chrUn_JTFH01001609v1_decoy": 1188, "chrUn_JTFH01001610v1_decoy": 1180, "chrUn_JTFH01001611v1_decoy": 1180, "chrUn_JTFH01001613v1_decoy": 1172, "chrUn_JTFH01001614v1_decoy": 1168, "chrUn_JTFH01001615v1_decoy": 1166, "chrUn_JTFH01001616v1_decoy": 1157, "chrUn_JTFH01001617v1_decoy": 1156, "chrUn_JTFH01001618v1_decoy": 1156, "chrUn_JTFH01001619v1_decoy": 1155, "chrUn_JTFH01001620v1_decoy": 1154, "chrUn_JTFH01001621v1_decoy": 1154, "chrUn_JTFH01001622v1_decoy": 1149, "chrUn_JTFH01001625v1_decoy": 1140, "chrUn_JTFH01001626v1_decoy": 1137, "chrUn_JTFH01001627v1_decoy": 1135, "chrUn_JTFH01001628v1_decoy": 1135, "chrUn_JTFH01001629v1_decoy": 1135, "chrUn_JTFH01001630v1_decoy": 1127, "chrUn_JTFH01001631v1_decoy": 1127, "chrUn_JTFH01001632v1_decoy": 1126, "chrUn_JTFH01001633v1_decoy": 1123, "chrUn_JTFH01001634v1_decoy": 1123, "chrUn_JTFH01001635v1_decoy": 1123, "chrUn_JTFH01001636v1_decoy": 1122, "chrUn_JTFH01001637v1_decoy": 1122, "chrUn_JTFH01001640v1_decoy": 1119, "chrUn_JTFH01001641v1_decoy": 1119, "chrUn_JTFH01001642v1_decoy": 1119, "chrUn_JTFH01001643v1_decoy": 1118, "chrUn_JTFH01001644v1_decoy": 1115, "chrUn_JTFH01001645v1_decoy": 1106, "chrUn_JTFH01001646v1_decoy": 1106, "chrUn_JTFH01001647v1_decoy": 1104, "chrUn_JTFH01001648v1_decoy": 1102, "chrUn_JTFH01001649v1_decoy": 1101, "chrUn_JTFH01001650v1_decoy": 1098, "chrUn_JTFH01001651v1_decoy": 1098, "chrUn_JTFH01001652v1_decoy": 1096, "chrUn_JTFH01001653v1_decoy": 1096, "chrUn_JTFH01001654v1_decoy": 1095, "chrUn_JTFH01001655v1_decoy": 1093, "chrUn_JTFH01001656v1_decoy": 1090, "chrUn_JTFH01001657v1_decoy": 1089, "chrUn_JTFH01001658v1_decoy": 1087, "chrUn_JTFH01001659v1_decoy": 1087, "chrUn_JTFH01001660v1_decoy": 1085, "chrUn_JTFH01001661v1_decoy": 1085, "chrUn_JTFH01001662v1_decoy": 1085, "chrUn_JTFH01001663v1_decoy": 1083, "chrUn_JTFH01001664v1_decoy": 1080, "chrUn_JTFH01001665v1_decoy": 1080, "chrUn_JTFH01001666v1_decoy": 1079, "chrUn_JTFH01001667v1_decoy": 1079, "chrUn_JTFH01001668v1_decoy": 1079, "chrUn_JTFH01001669v1_decoy": 1075, "chrUn_JTFH01001670v1_decoy": 1074, "chrUn_JTFH01001671v1_decoy": 1073, "chrUn_JTFH01001672v1_decoy": 1070, "chrUn_JTFH01001673v1_decoy": 1068, "chrUn_JTFH01001674v1_decoy": 1067, "chrUn_JTFH01001675v1_decoy": 1066, "chrUn_JTFH01001676v1_decoy": 1066, "chrUn_JTFH01001677v1_decoy": 1066, "chrUn_JTFH01001678v1_decoy": 1063, "chrUn_JTFH01001679v1_decoy": 1063, "chrUn_JTFH01001680v1_decoy": 1063, "chrUn_JTFH01001681v1_decoy": 1062, "chrUn_JTFH01001682v1_decoy": 1058, "chrUn_JTFH01001683v1_decoy": 1056, "chrUn_JTFH01001684v1_decoy": 1052, "chrUn_JTFH01001685v1_decoy": 1051, "chrUn_JTFH01001686v1_decoy": 1051, "chrUn_JTFH01001687v1_decoy": 1050, "chrUn_JTFH01001689v1_decoy": 1046, "chrUn_JTFH01001690v1_decoy": 1046, "chrUn_JTFH01001692v1_decoy": 1043, "chrUn_JTFH01001693v1_decoy": 1038, "chrUn_JTFH01001694v1_decoy": 1036, "chrUn_JTFH01001695v1_decoy": 1035, "chrUn_JTFH01001696v1_decoy": 1035, "chrUn_JTFH01001697v1_decoy": 1035, "chrUn_JTFH01001698v1_decoy": 1033, "chrUn_JTFH01001699v1_decoy": 1032, "chrUn_JTFH01001700v1_decoy": 1031, "chrUn_JTFH01001704v1_decoy": 1023, "chrUn_JTFH01001705v1_decoy": 1022, "chrUn_JTFH01001706v1_decoy": 1020, "chrUn_JTFH01001707v1_decoy": 1020, "chrUn_JTFH01001708v1_decoy": 1020, "chrUn_JTFH01001709v1_decoy": 1019, "chrUn_JTFH01001710v1_decoy": 1018, "chrUn_JTFH01001711v1_decoy": 1018, "chrUn_JTFH01001712v1_decoy": 1017, "chrUn_JTFH01001713v1_decoy": 1015, "chrUn_JTFH01001714v1_decoy": 1015, "chrUn_JTFH01001715v1_decoy": 1015, "chrUn_JTFH01001716v1_decoy": 1014, "chrUn_JTFH01001717v1_decoy": 1014, "chrUn_JTFH01001718v1_decoy": 1013, "chrUn_JTFH01001719v1_decoy": 1013, "chrUn_JTFH01001720v1_decoy": 1013, "chrUn_JTFH01001721v1_decoy": 1012, "chrUn_JTFH01001722v1_decoy": 1011, "chrUn_JTFH01001723v1_decoy": 1011, "chrUn_JTFH01001724v1_decoy": 1009, "chrUn_JTFH01001725v1_decoy": 1008, "chrUn_JTFH01001726v1_decoy": 1008, "chrUn_JTFH01001727v1_decoy": 1007, "chrUn_JTFH01001728v1_decoy": 1007, "chrUn_JTFH01001729v1_decoy": 1007, "chrUn_JTFH01001730v1_decoy": 1006, "chrUn_JTFH01001731v1_decoy": 1005, "chrUn_JTFH01001732v1_decoy": 1003, "chrUn_JTFH01001733v1_decoy": 1001, "chrUn_JTFH01001734v1_decoy": 1000, "chrUn_JTFH01001735v1_decoy": 19311, "chrUn_JTFH01001736v1_decoy": 11713, "chrUn_JTFH01001737v1_decoy": 11263, "chrUn_JTFH01001738v1_decoy": 9779, "chrUn_JTFH01001739v1_decoy": 9568, "chrUn_JTFH01001740v1_decoy": 9344, "chrUn_JTFH01001741v1_decoy": 9188, "chrUn_JTFH01001742v1_decoy": 9100, "chrUn_JTFH01001743v1_decoy": 8771, "chrUn_JTFH01001744v1_decoy": 8690, "chrUn_JTFH01001745v1_decoy": 8566, "chrUn_JTFH01001746v1_decoy": 8058, "chrUn_JTFH01001747v1_decoy": 7759, "chrUn_JTFH01001748v1_decoy": 7585, "chrUn_JTFH01001749v1_decoy": 7471, "chrUn_JTFH01001750v1_decoy": 7461, "chrUn_JTFH01001751v1_decoy": 7342, "chrUn_JTFH01001752v1_decoy": 7223, "chrUn_JTFH01001753v1_decoy": 7064, "chrUn_JTFH01001754v1_decoy": 6916, "chrUn_JTFH01001755v1_decoy": 6897, "chrUn_JTFH01001756v1_decoy": 6880, "chrUn_JTFH01001757v1_decoy": 6857, "chrUn_JTFH01001758v1_decoy": 6840, "chrUn_JTFH01001759v1_decoy": 6728, "chrUn_JTFH01001760v1_decoy": 6688, "chrUn_JTFH01001761v1_decoy": 6553, "chrUn_JTFH01001762v1_decoy": 6396, "chrUn_JTFH01001763v1_decoy": 6345, "chrUn_JTFH01001764v1_decoy": 6295, "chrUn_JTFH01001765v1_decoy": 6266, "chrUn_JTFH01001766v1_decoy": 6173, "chrUn_JTFH01001767v1_decoy": 6171, "chrUn_JTFH01001768v1_decoy": 6120, "chrUn_JTFH01001769v1_decoy": 6105, "chrUn_JTFH01001770v1_decoy": 6099, "chrUn_JTFH01001771v1_decoy": 5893, "chrUn_JTFH01001772v1_decoy": 5829, "chrUn_JTFH01001773v1_decoy": 5793, "chrUn_JTFH01001774v1_decoy": 5776, "chrUn_JTFH01001775v1_decoy": 5759, "chrUn_JTFH01001776v1_decoy": 5716, "chrUn_JTFH01001777v1_decoy": 5708, "chrUn_JTFH01001778v1_decoy": 5590, "chrUn_JTFH01001779v1_decoy": 5566, "chrUn_JTFH01001780v1_decoy": 5558, "chrUn_JTFH01001781v1_decoy": 5418, "chrUn_JTFH01001782v1_decoy": 5375, "chrUn_JTFH01001783v1_decoy": 5300, "chrUn_JTFH01001784v1_decoy": 5255, "chrUn_JTFH01001785v1_decoy": 5157, "chrUn_JTFH01001786v1_decoy": 5130, "chrUn_JTFH01001787v1_decoy": 4978, "chrUn_JTFH01001788v1_decoy": 4957, "chrUn_JTFH01001789v1_decoy": 4947, "chrUn_JTFH01001790v1_decoy": 4897, "chrUn_JTFH01001791v1_decoy": 4867, "chrUn_JTFH01001792v1_decoy": 4845, "chrUn_JTFH01001793v1_decoy": 4678, "chrUn_JTFH01001794v1_decoy": 4641, "chrUn_JTFH01001795v1_decoy": 4592, "chrUn_JTFH01001796v1_decoy": 4543, "chrUn_JTFH01001797v1_decoy": 4532, "chrUn_JTFH01001798v1_decoy": 4503, "chrUn_JTFH01001799v1_decoy": 4495, "chrUn_JTFH01001800v1_decoy": 4444, "chrUn_JTFH01001801v1_decoy": 4414, "chrUn_JTFH01001802v1_decoy": 4409, "chrUn_JTFH01001803v1_decoy": 4302, "chrUn_JTFH01001804v1_decoy": 4300, "chrUn_JTFH01001805v1_decoy": 4277, "chrUn_JTFH01001806v1_decoy": 4173, "chrUn_JTFH01001807v1_decoy": 4169, "chrUn_JTFH01001808v1_decoy": 4136, "chrUn_JTFH01001809v1_decoy": 4101, "chrUn_JTFH01001810v1_decoy": 4089, "chrUn_JTFH01001811v1_decoy": 4015, "chrUn_JTFH01001812v1_decoy": 4000, "chrUn_JTFH01001813v1_decoy": 3973, "chrUn_JTFH01001814v1_decoy": 3732, "chrUn_JTFH01001815v1_decoy": 3709, "chrUn_JTFH01001816v1_decoy": 3686, "chrUn_JTFH01001817v1_decoy": 3676, "chrUn_JTFH01001818v1_decoy": 3673, "chrUn_JTFH01001819v1_decoy": 3672, "chrUn_JTFH01001820v1_decoy": 3633, "chrUn_JTFH01001821v1_decoy": 3633, "chrUn_JTFH01001822v1_decoy": 3613, "chrUn_JTFH01001823v1_decoy": 3605, "chrUn_JTFH01001824v1_decoy": 3592, "chrUn_JTFH01001825v1_decoy": 3586, "chrUn_JTFH01001826v1_decoy": 3584, "chrUn_JTFH01001827v1_decoy": 3577, "chrUn_JTFH01001828v1_decoy": 3537, "chrUn_JTFH01001829v1_decoy": 3510, "chrUn_JTFH01001830v1_decoy": 3509, "chrUn_JTFH01001831v1_decoy": 3488, "chrUn_JTFH01001832v1_decoy": 3473, "chrUn_JTFH01001833v1_decoy": 3445, "chrUn_JTFH01001834v1_decoy": 3427, "chrUn_JTFH01001835v1_decoy": 3395, "chrUn_JTFH01001836v1_decoy": 3367, "chrUn_JTFH01001837v1_decoy": 3337, "chrUn_JTFH01001838v1_decoy": 3324, "chrUn_JTFH01001839v1_decoy": 3315, "chrUn_JTFH01001840v1_decoy": 3313, "chrUn_JTFH01001841v1_decoy": 3283, "chrUn_JTFH01001842v1_decoy": 3250, "chrUn_JTFH01001843v1_decoy": 3247, "chrUn_JTFH01001844v1_decoy": 3237, "chrUn_JTFH01001845v1_decoy": 3235, "chrUn_JTFH01001846v1_decoy": 3200, "chrUn_JTFH01001847v1_decoy": 3195, "chrUn_JTFH01001848v1_decoy": 3175, "chrUn_JTFH01001849v1_decoy": 3158, "chrUn_JTFH01001850v1_decoy": 3143, "chrUn_JTFH01001851v1_decoy": 3139, "chrUn_JTFH01001852v1_decoy": 3138, "chrUn_JTFH01001853v1_decoy": 3136, "chrUn_JTFH01001854v1_decoy": 3132, "chrUn_JTFH01001855v1_decoy": 3132, "chrUn_JTFH01001856v1_decoy": 3095, "chrUn_JTFH01001857v1_decoy": 3094, "chrUn_JTFH01001858v1_decoy": 3093, "chrUn_JTFH01001859v1_decoy": 3059, "chrUn_JTFH01001860v1_decoy": 2985, "chrUn_JTFH01001861v1_decoy": 2975, "chrUn_JTFH01001862v1_decoy": 2967, "chrUn_JTFH01001863v1_decoy": 2961, "chrUn_JTFH01001864v1_decoy": 2955, "chrUn_JTFH01001865v1_decoy": 2935, "chrUn_JTFH01001866v1_decoy": 2933, "chrUn_JTFH01001867v1_decoy": 2909, "chrUn_JTFH01001868v1_decoy": 2904, "chrUn_JTFH01001869v1_decoy": 2892, "chrUn_JTFH01001870v1_decoy": 2886, "chrUn_JTFH01001871v1_decoy": 2885, "chrUn_JTFH01001872v1_decoy": 2878, "chrUn_JTFH01001873v1_decoy": 2875, "chrUn_JTFH01001874v1_decoy": 2861, "chrUn_JTFH01001875v1_decoy": 2856, "chrUn_JTFH01001876v1_decoy": 2838, "chrUn_JTFH01001877v1_decoy": 2801, "chrUn_JTFH01001878v1_decoy": 2797, "chrUn_JTFH01001879v1_decoy": 2788, "chrUn_JTFH01001880v1_decoy": 2773, "chrUn_JTFH01001881v1_decoy": 2755, "chrUn_JTFH01001882v1_decoy": 2754, "chrUn_JTFH01001883v1_decoy": 2743, "chrUn_JTFH01001884v1_decoy": 2725, "chrUn_JTFH01001885v1_decoy": 2722, "chrUn_JTFH01001886v1_decoy": 2682, "chrUn_JTFH01001887v1_decoy": 2669, "chrUn_JTFH01001888v1_decoy": 2663, "chrUn_JTFH01001889v1_decoy": 2652, "chrUn_JTFH01001890v1_decoy": 2647, "chrUn_JTFH01001891v1_decoy": 2635, "chrUn_JTFH01001892v1_decoy": 2633, "chrUn_JTFH01001893v1_decoy": 2629, "chrUn_JTFH01001894v1_decoy": 2612, "chrUn_JTFH01001895v1_decoy": 2599, "chrUn_JTFH01001896v1_decoy": 2566, "chrUn_JTFH01001897v1_decoy": 2556, "chrUn_JTFH01001898v1_decoy": 2551, "chrUn_JTFH01001899v1_decoy": 2551, "chrUn_JTFH01001900v1_decoy": 2538, "chrUn_JTFH01001901v1_decoy": 2538, "chrUn_JTFH01001902v1_decoy": 2525, "chrUn_JTFH01001903v1_decoy": 2498, "chrUn_JTFH01001904v1_decoy": 2496, "chrUn_JTFH01001905v1_decoy": 2483, "chrUn_JTFH01001906v1_decoy": 2475, "chrUn_JTFH01001907v1_decoy": 2469, "chrUn_JTFH01001908v1_decoy": 2455, "chrUn_JTFH01001909v1_decoy": 2444, "chrUn_JTFH01001910v1_decoy": 2437, "chrUn_JTFH01001911v1_decoy": 2435, "chrUn_JTFH01001912v1_decoy": 2427, "chrUn_JTFH01001913v1_decoy": 2419, "chrUn_JTFH01001914v1_decoy": 2413, "chrUn_JTFH01001915v1_decoy": 2412, "chrUn_JTFH01001916v1_decoy": 2400, "chrUn_JTFH01001917v1_decoy": 2399, "chrUn_JTFH01001918v1_decoy": 2396, "chrUn_JTFH01001919v1_decoy": 2393, "chrUn_JTFH01001920v1_decoy": 2386, "chrUn_JTFH01001921v1_decoy": 2384, "chrUn_JTFH01001922v1_decoy": 2382, "chrUn_JTFH01001923v1_decoy": 2382, "chrUn_JTFH01001924v1_decoy": 2367, "chrUn_JTFH01001925v1_decoy": 2366, "chrUn_JTFH01001926v1_decoy": 2362, "chrUn_JTFH01001927v1_decoy": 2361, "chrUn_JTFH01001928v1_decoy": 2353, "chrUn_JTFH01001929v1_decoy": 2349, "chrUn_JTFH01001930v1_decoy": 2348, "chrUn_JTFH01001931v1_decoy": 2340, "chrUn_JTFH01001932v1_decoy": 2339, "chrUn_JTFH01001933v1_decoy": 2336, "chrUn_JTFH01001934v1_decoy": 2333, "chrUn_JTFH01001935v1_decoy": 2330, "chrUn_JTFH01001936v1_decoy": 2327, "chrUn_JTFH01001938v1_decoy": 2293, "chrUn_JTFH01001939v1_decoy": 2292, "chrUn_JTFH01001940v1_decoy": 2287, "chrUn_JTFH01001943v1_decoy": 2267, "chrUn_JTFH01001944v1_decoy": 2260, "chrUn_JTFH01001945v1_decoy": 2257, "chrUn_JTFH01001946v1_decoy": 2240, "chrUn_JTFH01001947v1_decoy": 2239, "chrUn_JTFH01001948v1_decoy": 2232, "chrUn_JTFH01001949v1_decoy": 2230, "chrUn_JTFH01001950v1_decoy": 2230, "chrUn_JTFH01001951v1_decoy": 2222, "chrUn_JTFH01001952v1_decoy": 2216, "chrUn_JTFH01001953v1_decoy": 2214, "chrUn_JTFH01001954v1_decoy": 2210, "chrUn_JTFH01001955v1_decoy": 2203, "chrUn_JTFH01001956v1_decoy": 2197, "chrUn_JTFH01001957v1_decoy": 2196, "chrUn_JTFH01001958v1_decoy": 2196, "chrUn_JTFH01001959v1_decoy": 2179, "chrUn_JTFH01001960v1_decoy": 2178, "chrUn_JTFH01001961v1_decoy": 2178, "chrUn_JTFH01001962v1_decoy": 2172, "chrUn_JTFH01001963v1_decoy": 2170, "chrUn_JTFH01001964v1_decoy": 2167, "chrUn_JTFH01001965v1_decoy": 2167, "chrUn_JTFH01001966v1_decoy": 2157, "chrUn_JTFH01001967v1_decoy": 2153, "chrUn_JTFH01001968v1_decoy": 2151, "chrUn_JTFH01001969v1_decoy": 2147, "chrUn_JTFH01001971v1_decoy": 2142, "chrUn_JTFH01001972v1_decoy": 2142, "chrUn_JTFH01001973v1_decoy": 2136, "chrUn_JTFH01001974v1_decoy": 2130, "chrUn_JTFH01001975v1_decoy": 2128, "chrUn_JTFH01001976v1_decoy": 2126, "chrUn_JTFH01001977v1_decoy": 2126, "chrUn_JTFH01001978v1_decoy": 2119, "chrUn_JTFH01001979v1_decoy": 2107, "chrUn_JTFH01001980v1_decoy": 2091, "chrUn_JTFH01001981v1_decoy": 2087, "chrUn_JTFH01001982v1_decoy": 2086, "chrUn_JTFH01001983v1_decoy": 2083, "chrUn_JTFH01001984v1_decoy": 2075, "chrUn_JTFH01001985v1_decoy": 2075, "chrUn_JTFH01001986v1_decoy": 2072, "chrUn_JTFH01001987v1_decoy": 2068, "chrUn_JTFH01001988v1_decoy": 2067, "chrUn_JTFH01001989v1_decoy": 2055, "chrUn_JTFH01001990v1_decoy": 2051, "chrUn_JTFH01001991v1_decoy": 2050, "chrUn_JTFH01001992v1_decoy": 2033, "chrUn_JTFH01001993v1_decoy": 2024, "chrUn_JTFH01001994v1_decoy": 2016, "chrUn_JTFH01001995v1_decoy": 2011, "chrUn_JTFH01001996v1_decoy": 2009, "chrUn_JTFH01001997v1_decoy": 2003, "chrUn_JTFH01001998v1_decoy": 2001}

# Mus musculus

mm7 = {
    "chr1": 194923535,
    "chr2": 182548267,
    "chrX": 164906252,
    "chr3": 159849039,
    "chr4": 155175443,
    "chr5": 153054177,
    "chr6": 149646834,
    "chr7": 141766352,
    "chr10": 130066766,
    "chr8": 127874053,
    "chr9": 123828236,
    "chr11": 122091587,
    "chr14": 119226840,
    "chr12": 117814103,
    "chr13": 116696528,
    "chr15": 103647385,
    "chr16": 98481019,
    "chr17": 93276925,
    "chr18": 90918714,
    "chr19": 61223509,
    "chrY": 15523453
}

mm8 = {
    "chr1": 197069962,
    "chr2": 181976762,
    "chr3": 159872112,
    "chr4": 155029701,
    "chr5": 152003063,
    "chr6": 149525685,
    "chr7": 145134094,
    "chr8": 132085098,
    "chr9": 124000669,
    "chrX": 165556469,
    "chrY": 16029404,
    "chr10": 129959148,
    "chr11": 121798632,
    "chr12": 120463159,
    "chr13": 120614378,
    "chr14": 123978870,
    "chr15": 103492577,
    "chr16": 98252459,
    "chr17": 95177420,
    "chr18": 90736837,
    "chr19": 61321190,
}


mm9 = {
    "chr1": 197195432,
    "chr2": 181748087,
    "chr3": 159599783,
    "chr4": 155630120,
    "chr5": 152537259,
    "chr6": 149517037,
    "chr7": 152524553,
    "chr8": 131738871,
    "chr9": 124076172,
    "chr10": 129993255,
    "chr11": 121843856,
    "chr12": 121257530,
    "chr13": 120284312,
    "chr14": 125194864,
    "chr15": 103494974,
    "chr16": 98319150,
    "chr17": 95272651,
    "chr18": 90772031,
    "chr19": 61342430,
    "chrX": 166650296,
    "chrY": 15902555
}

mm10 = {
    "chr1": 195471971,
    "chr2": 182113224,
    "chrX": 171031299,
    "chr3": 160039680,
    "chr4": 156508116,
    "chr5": 151834684,
    "chr6": 149736546,
    "chr7": 145441459,
    "chr10": 130694993,
    "chr8": 129401213,
    "chr14": 124902244,
    "chr9": 124595110,
    "chr11": 122082543,
    "chr13": 120421639,
    "chr12": 120129022,
    "chr15": 104043685,
    "chr16": 98207768,
    "chr17": 94987271,
    "chrY": 91744698,
    "chr18": 90702639,
    "chr19": 61431566
}

mm39 = {
    "chr1": 195154279,
    "chr2": 181755017,
    "chrX": 169476592,
    "chr3": 159745316,
    "chr4": 156860686,
    "chr5": 151758149,
    "chr6": 149588044,
    "chr7": 144995196,
    "chr10": 130530862,
    "chr8": 130127694,
    "chr14": 125139656,
    "chr9": 124359700,
    "chr11": 121973369,
    "chr13": 120883175,
    "chr12": 120092757,
    "chr15": 104073951,
    "chr16": 98008968,
    "chr17": 95294699,
    "chrY": 91455967,
    "chr18": 90720763,
    "chr19": 61420004,
}


# Drosophila melanogaster

dm5 = {
    "X": 22422827,
    "2L": 23011544,
    "2R": 21146708,
    "3L": 24543557,
    "3R": 27905053,
    "4": 1351857,
    "XHet": 204112,
    "2LHet": 368872,
    "2RHet": 3288761,
    "3LHet": 2555491,
    "3RHet": 2517507,
    "YHet": 347038,
    "MT": 19517
}

dm6 = {
    "X": 23542271,
    "2L": 23513712,
    "2R": 25286936,
    "3L": 28110227,
    "3R": 32079331,
    "4": 1348131,
    "Y": 3667352
}

# Danio Rerio

danRer10 = {
    "Chromosome 1": 58871917,
    "Chromosome 2": 59543403,
    "Chromosome 3": 62385949,
    "Chromosome 4": 76625712,
    "Chromosome 5": 71715914,
    "Chromosome 6": 60272633,
    "Chromosome 7": 74082188,
    "Chromosome 8": 54191831,
    "Chromosome 9": 56892771,
    "Chromosome 10": 45574255,
    "Chromosome 11": 45107271,
    "Chromosome 12": 49229541,
    "Chromosome 13": 51780250,
    "Chromosome 14": 51944548,
    "Chromosome 15": 47771147,
    "Chromosome 16": 55381981,
    "Chromosome 17": 53345113,
    "Chromosome 18": 51008593,
    "Chromosome 19": 48790377,
    "Chromosome 20": 55370968,
    "Chromosome 21": 45895719,
    "Chromosome 22": 39226288,
    "Chromosome 23": 46272358,
    "Chromosome 24": 42251103,
    "Chromosome 25": 36898761
}

danRer11 = {
    "Chromosome 1": 59578282,
    "Chromosome 2": 59640629,
    "Chromosome 3": 62628489,
    "Chromosome 4": 78093715,
    "Chromosome 5": 72500376,
    "Chromosome 6": 60270059,
    "Chromosome 7": 74282399,
    "Chromosome 8": 54304671,
    "Chromosome 9": 56459846,
    "Chromosome 10": 45420867,
    "Chromosome 11": 45484837,
    "Chromosome 12": 49182954,
    "Chromosome 13": 52186027,
    "Chromosome 14": 52660232,
    "Chromosome 15": 48040578,
    "Chromosome 16": 55266484,
    "Chromosome 17": 53461100,
    "Chromosome 18": 51023478,
    "Chromosome 19": 48449771,
    "Chromosome 20": 55201332,
    "Chromosome 21": 45934066,
    "Chromosome 22": 39133080,
    "Chromosome 23": 46223584,
    "Chromosome 24": 42172926,
    "Chromosome 25": 37502051
}

# Caenorhabditis elegans

WBcel215 = {
    "Chromosome I": 15072423,
    "Chromosome II": 15279345,
    "Chromosome III": 13783700,
    "Chromosome IV": 17493793,
    "Chromosome V": 20924149,
    "Chromosome X": 17718866
}

WBcel235 = {
    "I": 15072434,
    "II": 15279421,
    "III": 13783801,
    "IV": 17493829,
    "V": 20924180,
    "X": 17718942,
    "MT": 13794
}

# Rattus norvegicus

mRatBN7_2 = {
    "1": 260522016,
    "2": 249053267,
    "3": 169034231,
    "4": 182687754,
    "5": 166875058,
    "6": 140994061,
    "7": 135012528,
    "8": 123900184,
    "9": 114175309,
    "10": 107211142,
    "11": 86241447,
    "12": 46669029,
    "13": 106807694,
    "14": 104886043,
    "15": 101769107,
    "16": 84729064,
    "17": 86533673,
    "18": 83828827,
    "19": 57337602,
    "20": 54435887,
    "X": 152453651,
    "Y": 18315841
}

Rnor_6_0 = {
    "1": 282763074,
    "2": 266435125,
    "3": 177699992,
    "4": 184226339,
    "5": 173707219,
    "6": 147991367,
    "7": 145729302,
    "8": 133307652,
    "9": 122095297,
    "10": 112626471,
    "11": 90463843,
    "12": 52716770,
    "13": 114033958,
    "14": 115493446,
    "15": 111246239,
    "16": 90668790,
    "17": 90843779,
    "18": 88201929,
    "19": 62275575,
    "20": 56205956,
    "X": 159970021,
    "Y": 3310458
}

# Saccharomyces cerevisiae

R64 = {
    "I": 230218,
    "II": 813184,
    "III": 316620,
    "IV": 1531933,
    "V": 576874,
    "VI": 270161,
    "VII": 1090940,
    "VIII": 562643,
    "IX": 439888,
    "X": 745751,
    "XI": 666816,
    "XII": 1078177,
    "XIII": 924431,
    "XIV": 784333,
    "XV": 1091291,
    "XVI": 948066,
    "MT": 85779
}

# E. coli

ASM584v2 = {
    "GCF_000005845.2": 4641652
}

ASM886v2 = {
    "chromosome": 2778819,
    "pOSAK1": 1435,
    "pO157": 44135
}

# Pan troglodytes

pantro3_0 = {
    "1": 228573443,
    "2A": 111504155,
    "2B": 133216015,
    "3": 202621043,
    "4": 194502333,
    "5": 181907262,
    "6": 175400573,
    "7": 166211670,
    "8": 147911612,
    "9": 116767853,
    "10": 135926727,
    "11": 135753878,
    "12": 137163284,
    "13": 100452976,
    "14": 91965084,
    "15": 83230942,
    "16": 81586097,
    "17": 83181570,
    "18": 78221452,
    "19": 61309027,
    "20": 66533130,
    "21": 33445071,
    "22": 37823149,
    "X": 155549662,
    "Y": 26350515
}

pantro_2_1_4 = {
    "1": 228333871,
    "2A": 113622374,
    "2B": 247518478,
    "3": 202329955,
    "4": 193495092,
    "5": 182651097,
    "6": 172623881,
    "7": 161824586,
    "8": 143986469,
    "9": 137840987,
    "10": 133524379,
    "11": 133121534,
    "12": 134246214,
    "13": 115123233,
    "14": 106544938,
    "15": 99548318,
    "16": 89983829,
    "17": 82630442,
    "18": 76611499,
    "19": 63644993,
    "20": 61729293,
    "21": 32799110,
    "22": 49737984,
    "X": 156848144,
    "Y": 23952694
}


# Macaca mulatta

Mmul10 = {
    "1": 223616942,
    "2": 196197964,
    "3": 185288947,
    "4": 169963040,
    "5": 187317192,
    "6": 179085566,
    "7": 169868564,
    "8": 145679320,
    "9": 134124166,
    "10": 99517758,
    "11": 133066086,
    "12": 130043856,
    "13": 108737130,
    "14": 128056306,
    "15": 113283604,
    "16": 79627064,
    "17": 95433459,
    "18": 74474043,
    "19": 58315233,
    "20": 77137495,
    "X": 153388924,
    "Y": 11753682
}

rheMac8 = {
    "chr1": 225584828,
    "chr2": 204787373,
    "chr5": 190429646,
    "chr3": 185818997,
    "chr6": 180051392,
    "chr4": 172585720,
    "chr7": 169600520,
    "chrX": 149150640,
    "chr8": 144306982,
    "chr11": 133663169,
    "chr9": 129882849,
    "chr14": 127894412,
    "chr12": 125506784,
    "chr15": 111343173,
    "chr13": 108979918,
    "chr17": 95684472,
    "chr10": 92844088,
    "chr16": 77216781,
    "chr20": 74971481,
    "chr18": 70235451,
    "chr19": 53671032,
    "chrY": 11753682,
}

rheMac3 = {
    "1": 229590362,
    "2": 192599291,
    "3": 198365852,
    "4": 169853907,
    "5": 184069350,
    "6": 180491593,
    "7": 170124641,
    "8": 150158102,
    "9": 130959383,
    "10": 95678458,
    "11": 134985053,
    "12": 106825570,
    "13": 137160470,
    "14": 135538710,
    "15": 112192347,
    "16": 80904711,
    "17": 94676912,
    "18": 74429130,
    "19": 64156825,
    "20": 88476065,
    "X": 155416530
}

# Arabidopsis thaliana

TAIR = {
    "GCF_000001735.4_1": 30427671,
    "GCF_000001735.4_2": 19698289,
    "GCF_000001735.4_3": 23459830,
    "GCF_000001735.4_4": 18585056,
    "GCF_000001735.4_5": 26975502,
    "GCF_000001735.4_6": 367808,
    "GCF_000001735.4_7": 154478
}

# Sus scrofa

Sscrofa10_2 = {
    "chr1": 315321322,
    "chr2": 162569375,
    "chr3": 144787322,
    "chr4": 143465943,
    "chr5": 111506441,
    "chr6": 157765593,
    "chr7": 134764511,
    "chr8": 148491826,
    "chr9": 153670197,
    "chr10": 79102373,
    "chr11": 87690581,
    "chr12": 63588571,
    "chr13": 218635234,
    "chr14": 153851969,
    "chr15": 157681621,
    "chr16": 86898991,
    "chr17": 69701581,
    "chr18": 61220071,
    "chrX": 144288218,
    "chrY": 1637716
}

Sscrofa11_1 = {
    "chr1": 274330532,
    "chr2": 151935994,
    "chr3": 132848913,
    "chr4": 130910915,
    "chr5": 104526007,
    "chr6": 170843587,
    "chr7": 121844099,
    "chr8": 138966237,
    "chr9": 139512083,
    "chr10": 69359453,
    "chr11": 79169978,
    "chr12": 61602749,
    "chr13": 208334590,
    "chr14": 141755446,
    "chr15": 140412725,
    "chr16": 79944280,
    "chr17": 63494081,
    "chr18": 55982971,
    "chrX": 125939595
}

## Lists and dictionaries for comparison()
#major_releases = {"hg16": hg16, "hg17": hg17, "hg18": hg18, "GRCh37": GRCh37, "GRCh38": GRCh38, "T2T": T2T}

major_releases = {"hg16": {"ref_gen": hg16, "build": "hg16", "species": "Homo sapiens"},
                  "hg17": {"ref_gen": hg17, "build": "hg17", "species": "Homo sapiens"},
                  "hg18": {"ref_gen": hg18, "build": "hg18", "species": "Homo sapiens"},
                  "GRCh37": {"ref_gen": GRCh37, "build": "GRCh37", "species": "Homo sapiens"},
                  "GRCh38": {"ref_gen": GRCh38, "build": "GRCh38", "species": "Homo sapiens"},
                  "T2T": {"ref_gen": T2T, "build": "T2T", "species": "Homo sapiens"},
                  "mm7": {"ref_gen": mm7, "build": "mm7", "species": "Mus Musculus"},
                  "mm8": {"ref_gen": mm8, "build": "mm8", "species": "Mus musculus"},
                  "mm9": {"ref_gen": mm9, "build": "mm9", "species": "Mus musculus"},
                  "mm10": {"ref_gen": mm10, "build": "mm10", "species": "Mus musculus"},
                  "mm39": {"ref_gen": mm39, "build": "mm39", "species": "Mus musculus"},
                  "dm5": {"ref_gen": dm5, "build": "dm5", "species": "Drosophila Melanogaster"},
                  "dm6": {"ref_gen": dm6, "build": "dm6", "species": "Drosophila Melanogaster"},
                  "danRer10": {"ref_gen": danRer10, "build": "danRer10", "species": "Danio Rerio"},
                  "danRer11": {"ref_gen": danRer11, "build": "danRer11", "species": "Danio Rerio"},
                  "WBcel215": {"ref_gen": WBcel215, "build": "WBcel215", "species": "Caenorhabditis elegans"},
                  "WBcel235": {"ref_gen": WBcel235, "build": "WBcel235", "species": "Caenorhabditis elegans"},
                  "mRatBN7_2": {"ref_gen": mRatBN7_2, "build": "mRatBN7_2", "species": "Rattus norvegicus"},
                  "Rnor_6_0": {"ref_gen": Rnor_6_0, "build": "Rnor_6_0", "species": "Rattus norvegicus"},
                  "R64": {"ref_gen": R64, "build": "R64", "species": "Saccharomyces cerevisiae"},
                  "ASM886v2": {"ref_gen": ASM886v2, "build": "ASM886v2", "species": "Escherichia coli"},
                  "ASM584v2": {"ref_gen": ASM584v2, "build": "ASM584v2", "species": "Escherichia Coli"},
                  "pantro3_0": {"ref_gen": pantro3_0, "build": "pantro3_0", "species": "Pan troglodytes"},
                  "pantro_2_1_4": {"ref_gen": pantro_2_1_4, "build": "Pan_troglodytes-2.1",
                                   "species": "Pan troglodytes"},
                  "Mmul10": {"ref_gen": Mmul10, "build": "Mmul10", "species": "Macaca mulatta"},
                  "rheMac8": {"ref_gen": rheMac8, "build": "rheMac8", "species": "Macaca mulatta"},
                  "rheMac3": {"ref_gen": rheMac3, "build": "rheMac3", "species": "Macaca mulatta"},
                  "TAIR": {"ref_gen": TAIR, "build": "TAIR", "species": "Arabidopsis thaliana"},
                  "Sscrofa10_2": {"ref_gen": Sscrofa10_2, "build": "Sscrofa10_2", "species": "Sus scrofa"},
                  "Sscrofa11_1": {"ref_gen": Sscrofa11_1, "build": "Sscrofa11_1", "species": "Sus scrofa"}
                  }


flavors_GRCh37 = {"b37":  {"ref_gen": b37, "build": "b37", "species": "Homo sapiens"},
                  "hs37d5": {"ref_gen": hs37d5, "build": "hs37d5", "species": "Homo sapiens"},
                  "hg19":  {"ref_gen": hg19, "build": "hg19", "species": "Homo sapiens"}
                  }  # GRCh37 flavors that can be inferred

"""
IMPORTANT NOTE
The order of the keys in flavors_GRCh37 is not modifiable. 
b37 must be before hs37d5 for when it's
"""

def get_duplicate_lengths():
    nested_list = []

    for key, value in major_releases.items():
        ref_gen_dict = value["ref_gen"]
        nested_list.append(list(ref_gen_dict.values()))  # Convert dict values to a list before appending

    values_list = []

    for i in nested_list:
        for j in i:
            values_list.append(j)

    from collections import Counter

    # Count occurrences
    counts = Counter(values_list)

    # Filter duplicates
    duplicates = {item: count for item, count in counts.items() if count > 1}

    print(duplicates)


