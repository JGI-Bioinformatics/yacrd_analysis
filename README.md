Summary:
In a sampling of 20 reads that yet another chimeric read detector (YACRD) detected in a run of E. coli (run X0161), 20 are true positives. The same results are seen on a metagenome run of Actino Mock (run X0163). No false positives are seen. The rate of true and false negatives is unknown. The github repository for YACRD can be viewed here: https://github.com/natir/yacrd. The repository for dictionaryChimeric.py is here: https://github.com/JGI-Bioinformatics/yacrd_analysis. 

Explanation of code:
The output of dictionaryChimeric.py is a new fastq file with all chimeric reads split out into separate reads. This new fastq file can now be used for analysis. The program takes the output of running YACRD sent to a txt file. This YACRD text file is parsed and stored in a dictionary where the read name is the key and the value is the rest of the output from YACRD for that read.
For each line in the fastq file that is to be used, a check is done to see if that line what flagged by YACRD as being chimeric or "not covered". If it is anything other than chimeric, nothing is done (i.e. no reads are written to the new fastq file that is being generated). If it is chimeric, then the read is split up into individual reads where there is sufficient coverage, and these reads are written to the new fastq file. If the line in the original fastq file was not flagged by YACRD, then the whole read is written to the new fastq file as is, because it is a valid read and nothing needs to be done.
Below this page is data from the E. coli run analysis of YACRD's impact on true/false positives, as well as that of Actino Mock.

Running the program:
	In order to run dictionaryChimeric.py, you will first need to run minimap2 and then YACRD. The commands I used are as follows: 
Run minimap2: 
shifter --image=robegan21/minimap2 minimap2 -x ava-ont X0163_small.fastq X0163_small.fastq > smallMapping.paf
X0163_small.fastq is the original fastq
smallMapping.paf is the mapping file I want to create
Run yacrd:
yacrd -f smallMapping.paf -i smallMapping.paf -o output.fastq > smallYacrd.txt
The paf file must be the paf file you created when running minimap2
Run dictionaryChimeric.py on NERSC:
source activate yacrd
conda install biopython
python dictionaryChimeric.py > smallOutput.txt
smallOutput.txt is a file that contains useful printed information, including each read's name, the number of pieces that the read will be separated into, the YACRD data for each piece (i.e. length, starting index, and ending index of each region), and the name of each new read
The file smallOutput.fastq is the resulting fastq file that can now be used for analysis
	
Results from E. coli run (X0161):
on newSmallFastq:
380,000 reads
4296 are chimeric
TP = true positive
FP = false positive

Chimeric read: 50a2b172-651f-4333-ae6d-7f2bda586a80_Basecall_1D_template
3400713
152669
ref distance: 3248044
read distance: 6253-6249 = 4
difference: 3248040
2 pieces
conclusion: TP
Chimeric read: e0206847-2acb-4d7b-9844-26d218d6aa9b_Basecall_1D_template
1111429
429948
ref distance: 681481
read distance: 385-378 = 7
difference: 681474
2 pieces
conclusion: TP
Chimeric read: 2e95da4b-9d51-4eaa-a1cb-815ae3724469_Basecall_1D_template
3888829
3889169
ref distance: 340
read distance: 700-694 = 6
difference: 334
2 pieces
conclusion: TP
Chimeric read: 4c3f327f-5ff8-4771-a9e6-9dfd4bfc2033_Basecall_1D_template
3243801
4068657
ref distance: 824856
read distance: 435-417 = 18
difference: 824838
2 pieces
conclusion: TP
Chimeric read: d8ec6051-ab6b-4cf6-9a6e-d4fb8dd1b05e_Basecall_1D_template
3400711
2407703
ref distance: 993008
read distance: 21577-20756 = 821
difference: 992187
2 pieces
conclusion: TP
Chimeric read: bb09585b-55b8-44ad-b929-e0e22e776035_Basecall_1D_template
504434
3312064
ref distance: 2807630
read distance: 514-385 = 129
difference: 2807501
2 pieces
conclusion: TP
Chimeric read: f6a24bd9-19a7-497c-8542-8834a903a9c2_Basecall_1D_template
1077570
1078140
ref distance: 570
read distance: 882-873 = 9
difference: 561
2 pieces
conclusion: FP?
Chimeric read: 328c98ff-cb8e-40bd-837a-d987f41ebbcf_Basecall_1D_template
3954929
2044402
ref distance: 1910527
read distance: 635-627 = 8
difference: 1910519
2 pieces
conclusion: TP
Chimeric read: 14c0a8fe-ec0f-4e6f-bb6b-f5466fa8f2be_Basecall_1D_template
4096222
1031137
ref distance: 3065085
read distance: 1561-1545 = 16
difference: 3065069
2 pieces
conclusion: TP
Chimeric read: 5285c16c-c99c-451e-93ce-ecb87d82255d_Basecall_1D_template
4614982
770017
ref distance: 3844965
read distance: 3160-3138 = 22
difference: 3844943
2 pieces
conclusion: TP
Chimeric read: c8008c2d-f2fe-4656-a42b-5affcf10b5fe_Basecall_1D_template
2975146
2321423
ref distance: 653723
read distance: 3960-3947 = 13
difference: 653710
2 pieces
conclusion: TP
Chimeric read: 727fcf64-11e7-4b65-b5ee-4cd351c2a0d2_Basecall_1D_template
2703702
3967537
ref distance: 1263835
read distance: 785-784 = 1
difference: 1263834
2 pieces
conclusion: TP
Chimeric read: b6584aec-5d4e-4b84-a973-d80f6b280d43_Basecall_1D_template
1407235
779166
ref distance: 628069
read distance: 541-539 = 2
difference: 628067
2 pieces
conclusion: TP
Chimeric read: edc14709-0652-4317-b017-0cfd6fab2494_Basecall_1D_template
3799432
3799667
ref distance: 235
read distance: 18610-18601 = 9
difference: 226
2 pieces
conclusion: FP?
Chimeric read: 8414a449-2caa-4ece-a06d-632e1306c326_Basecall_1D_template
2391354
1614167
ref distance: 777187
read distance: 66316-66277 = 39
difference: 777148
2 pieces
conclusion: TP
Chimeric read: 2b4fea49-75e0-4911-bcdb-2cb7507fa1b5_Basecall_1D_template
177889
178525
ref distance: 636
read distance: 3766-3761 = 5
difference: 631
2 pieces
conclusion: TP
Chimeric read: bf10b8f5-276f-4912-b33c-b9d79096a8dd_Basecall_1D_template
842171
833019
ref distance: 9152
read distance: 11589-11233 = 356
difference: 8796
2 pieces
conclusion: TP
Chimeric read: a0c7035c-ee8a-42d6-a65c-587ea96031bd_Basecall_1D_template
4351214
4351239
ref distance: 25
read distance: 1633-1632 = 1
difference: 24
2 pieces
conclusion: TP
Chimeric read: fabd41fd-e499-441e-84f8-a0924740cfe7_Basecall_1D_template
3598627
3414887
ref distance: 183740
read distance: 1354-1353 = 1
difference: 183739
2 pieces
conclusion: TP
Chimeric read: 0afb4f6f-d60d-4099-80b8-668056d4f9d2_Basecall_1D_template
756616
3114568
ref distance: 2357970
read distance: 176-171 = 5
difference: 2357965
2 pieces
conclusion: TP

Results from metagenome run (ActinoMock X0163):
/global/projectb/scratch/ashleigh/chimeric/X0163/X0163_small.fastq
500,000 reads
10,337 are chimeric

Chimeric read: b5039ed2-a9b3-4d9c-a435-c58a95055ea2_Basecall_1D_template
Position 1: 165022
Position 2: 3025296
Reference distance: 2860274
Read distance: 3539-3510 = 29
Difference: 2860245
Number of pieces: 2
Conclusion: TP
Chimeric read: bbcc16f5-5257-4dd3-9632-f6a7ed614886_Basecall_1D_template
Position 1: 2877556
Position 2: 255
Reference distance: 2877301
Read distance: 1282-1251 = 31
Difference: 2877270
Number of pieces: 2
Conclusion: TP
Chimeric read: 093571b1-d3a2-45d6-8ed2-db69ab472ced_Basecall_1D_template
Position 1: 2643725
Position 2: 171207
Reference distance: 2472518
Read distance: 885-824 = 61
Difference: 2472457
Number of pieces: 2
Conclusion: TP
Chimeric read: d8b2fb8e-0476-4f64-af71-a726a9ac9011_Basecall_1D_template
Position 1: 1044000
Position 2: 255
Reference distance: 1043745
Read distance: 1337-1311 = 26
Difference: 1043719
Number of pieces: 2
Conclusion: TP
Chimeric read: 61437a9a-e2b9-417e-99e0-02c809711586_Basecall_1D_template
Position 1: 625910
Position 2: 1149771
Reference distance: 523861
Read distance: 1705-1655 = 50
Difference: 523811
Number of pieces: 2
Conclusion: TP
Chimeric read: 11886edb-a844-47df-8f5e-8519d8d32ffb_Basecall_1D_template
Position 0: 1844594 (AM 8)
Position 1: 3386543 (AM 3)
Position 2: 1327974 (AM3)
Reference distance:
Read distance:1502-1499 = 3
Difference:
Number of pieces: 3
Conclusion: TP, definitely chimeric
Chimeric read: 7e7a9ba6-09f8-49c8-aa61-92410595f534_Basecall_1D_template
Position 1: 3990998
Position 2: 1816449
Reference distance: 2174549
Read distance:3532-3469=63
Difference: 2174486
Number of pieces:2
Conclusion: TP
Chimeric read: 1c423dd6-5db1-4ff5-879b-f37a7c9e0a37_Basecall_1D_template
Position 0: 537113 AM1a
Position 1: 254262 AM1a
Position 2: 229153 AM2
Reference distance:
Read distance:
Difference:
Number of pieces: 3
Conclusion: TP
Chimeric read: 92468353-8707-4c6d-babd-2e659da427ca_Basecall_1D_template
Position 1: 2655468
Position 2: 3456517
Reference distance: 801049
Read distance: 3289-3287 = 2
Difference: 801047
Number of pieces: 2
Conclusion: TP
Chimeric read: c537fffe-cd6d-4558-8421-2e2247364ccc_Basecall_1D_template
Position 1: 336712
Position 2: 332822
Reference distance: 3890
Read distance: 2674-2646 = 28
Difference: 3862
Number of pieces: 2
Conclusion: TP
Chimeric read: e014d7b6-8178-4d36-9267-39b6edf68e40_Basecall_1D_template
Position 1: 1067315
Position 2: 5289438
Reference distance: 4222123
Read distance: 799-744 = 55
Difference: 4222068
Number of pieces: 2
Conclusion: TP
Chimeric read: 069bae3d-0ec3-40f8-83aa-06d02639db7d_Basecall_1D_template
Position 1: 879318
Position 2: 881974
Reference distance: 2656
Read distance:4362-2513 = 1849
Difference: 807
Number of pieces: 2
Conclusion: TP
Chimeric read: 9643e482-426d-4c38-a3da-f77a8dd64a9f_Basecall_1D_template
Position 1: 2391165
Position 2: 1890673
Reference distance: 500492
Read distance: 3015-2990 = 25
Difference: 500467
Number of pieces: 2
Conclusion: TP
Chimeric read: a04f87cc-732d-4d70-9214-4db7c35aafda_Basecall_1D_template
Position 1: 1392146
Position 2: 1392139
Position 3: 3649554
Reference distance: 2257415 ? (3-1)
Read distance: 5023-1918 = 3105
Difference: 2254310
Number of pieces: 2 with small gap
Conclusion: TP, probably just a gap
Chimeric read: 50ef1361-a8be-40bd-8fc6-d05fc7f88ada_Basecall_1D_template
Position 1: 2744519
Position 2: 1020564
Reference distance: 1723955
Read distance:12712-12592 = 120
Difference: 1723835
Number of pieces: 2
Conclusion: TP
Chimeric read: 3406624d-aef2-49ca-8f50-2c23bbed7032_Basecall_1D_template
Position 1: 707114
Position 2: 2069936
Reference distance: 1362822
Read distance:3469-3415 = 54
Difference: 1362768
Number of pieces: 2
Conclusion: TP
Chimeric read: b7e8c4b5-1fbc-491a-b491-de1419f65b58_Basecall_1D_template
Position 1: 2516837
Position 2: 1127938
Reference distance: 1388899
Read distance:894-835 = 59
Difference: 1388840
Number of pieces:2
Conclusion: TP
Chimeric read: 04cef3b5-581b-41b1-b744-bd363f426747_Basecall_1D_template
Position 1: 770492
Position 2: 2031873
Reference distance: 1261381
Read distance:1517-1461=56
Difference: 1261325
Number of pieces:2
Conclusion: TP
Chimeric read: 97aaf4df-23de-4f4f-bd6a-91451b24491a_Basecall_1D_template
Position 1: 348181
Position 2: 1177910
Reference distance: 829729
Read distance:2125-2121=4
Difference: 829725
Number of pieces:2
Conclusion: TP
Chimeric read: 1225c4bf-7c20-49e8-9695-d0d75a3ceec4_Basecall_1D_template
Position 1: 2227483
Position 2: 1552265
Reference distance: 675218
Read distance:1785-1708=77
Difference: 675141
Number of pieces:2
Conclusion: TP
Analysis of Not_covered in X0163 (metagenome)
c5f22007-690c-4ecc-bf7e-e03acba85c59_Basecall_1D_template
Reference mapping: 0
a392a9d7-d46a-465f-895a-e072d6b7222d_Basecall_1D_template
Reference mapping: 0
2743687f-6256-4e93-a81a-fc2896df76f8_Basecall_1D_template:
Reference mapping: 0
4ef27459-592c-4a78-a5db-1c58205eb841_Basecall_1D_template
Reference mapping: 0
1c82aedb-9913-493b-a4b2-11a9ff19e158_Basecall_1D_template
Reference mapping: 0
d89a2cc2-9009-4de6-85a4-5a10d5502c88_Basecall_1D_template
Reference mapping: 0
b3e41bf7-3dce-4be1-a122-5002e4bf4a3b_Basecall_1D_template
Reference mapping: 0
89809a7a-95d5-4434-b6ab-9f50361138e2_Basecall_1D_template
Reference mapping: 0
4b5ba5ed-68c8-44e4-963d-8ec4e2c85928_Basecall_1D_template
Reference mapping: 0
8c587974-2d8e-4ad8-82a4-08c326833ded_Basecall_1D_template
Reference mapping: 0
2005fc05-4e15-4b40-9d39-3847b012e012_Basecall_1D_template
Reference mapping: 0
a7eb81ba-2f89-44cd-afd3-afa1a03dd716_Basecall_1D_template
Reference mapping: 0
b29f43f6-7c96-45d2-81ea-26e0b18d69e5_Basecall_1D_template
Reference mapping: 4392406
88b6e8dd-5845-4702-9bff-1428784f0131_Basecall_1D_template
Reference mapping: 0
90ef8a26-419c-4cf0-9fb5-169c8a0c54c2_Basecall_1D_template
Reference mapping: 0
e59ac1c6-e0f9-485b-b8b0-7fd85ba30b87_Basecall_1D_template
Reference mapping: 0
47103d13-46a9-4f27-afe5-d88c521c1414_Basecall_1D_template
Reference mapping: 0
f6dacf34-00c3-433c-bbc0-d031a18858a1_Basecall_1D_template
Reference mapping: 0
9ad90ff1-dd36-4e51-b7e9-57490f8694b0_Basecall_1D_template
Reference mapping: 0


