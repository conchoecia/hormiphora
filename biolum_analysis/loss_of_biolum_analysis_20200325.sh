#!/bin/bash

cd ~/genomes/hormiphora_californensis/hormiphora/annotation/Hcv1a1d20200325_release
#gzip -dc Hcv1a1d20200325.gff.gz > Hcv1a1d20200325.gff
#gzip -dc Hcv1a1d20200325_model_proteins.pep.gz > Hcv1a1d20200325_model_proteins.pep
gzip -dc Hcv1a1d20200325_transcripts.fasta.gz > Hcv1a1d20200325_transcripts.fasta
~/diamond-latest/diamond makedb --in Hcv1a1d20200325_model_proteins.pep --db Hcv1a1d20200325_model_proteins.pep

cd ~/genomes/hormiphora_californensis/hormiphora/biolum_analysis/
ln -s ../annotation/Hcv1a1d20200325_release/Hcv1a1d20200325_model_proteins.pep
ln -s ../annotation/Hcv1a1d20200325_release/Hcv1a1d20200325_model_proteins.pep.dmnd

ln -s ../synteny/Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab
ln -s ../synteny/Hcv1a1d20200325_model_proteins.pep.vs_ml2unfiltered.tab


# for the FYY genes
echo "# FYY searches" >> Hcv1a1d20200325_biolum_results.txt
grep ML199826a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c6.g634.i1	ML199826a	33.0	324	194	10	8	320	35	346	3.0e-48	189.1
#Hcv1.1.c11.g551.i1	ML199826a	42.3	324	171	4	47	356	34	355	8.3e-71	264.2
echo "# FYY neighboring genes" >> Hcv1a1d20200325_biolum_results.txt
grep ML199825a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c2.g480.i1	ML199825a	58.2	361	75	5	1	349	119	415	3.2e-104	376.7
grep ML199827a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c2.g480.i1	ML199827a	74.6	465	118	0	573	1037	1	465	9.0e-208	720.7
#Hcv1.1.c5.g538.i1	ML199827a	26.3	95	64	2	271	360	312	405	2.7e-04	43.5
grep ML199828a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c2.g481.i2	ML199828a	35.6	317	163	9	1	301	44	335	1.2e-35	147.1
#Hcv1.1.c9.g512.i1	ML199828a	51.0	49	24	0	139	187	202	250	9.4e-06	47.4


grep ML263511a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c6.g219.i1	ML263511a	36.3	543	327	8	37	568	71	605	1.6e-81	300.4
grep ML263512a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c6.g218.i1	ML263512a	31.0	542	361	6	46	578	89	626	1.1e-69	261.2
#Hcv1.1.c2.g778.i2	ML263512a	22.1	353	232	13	5	346	18	338	1.5e-04	45.1
grep ML263513a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c6.g217.i1	ML263513a	66.3	276	91	1	37	310	30	305	6.0e-112	401.4


grep ML263515a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c2.g345.i1	ML263515a	65.7	353	76	3	87	410	1	337	1.5e-129	459.5
#Hcv1.1.c6.g622.i1	ML263515a	55.4	352	115	3	80	404	1	337	1.7e-106	382.9
#Hcv1.1.c2.g700.i2	ML263515a	25.6	312	199	3	109	397	20	321	7.2e-28	121.7
#Hcv1.1.c3.g629.i1	ML263515a	28.8	330	215	6	80	394	2	326	1.2e-27	120.9
#Hcv1.1.c6.g172.i1	ML263515a	24.4	308	214	4	131	426	22	322	8.1e-20	95.5
#Hcv1.1.c12.g332.i1	ML263515a	27.2	287	196	5	126	399	42	328	3.0e-18	89.7
grep ML263516a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c2.g345.i1	ML263516a	83.5	139	23	0	272	410	1	139	3.9e-66	248.8
#Hcv1.1.c6.g622.i1	ML263516a	70.3	138	41	0	267	404	2	139	9.6e-57	217.6

echo "# against unfiltered set" >> Hcv1a1d20200325_biolum_results.txt
grep MLRB263543 Hcv1a1d20200325_model_proteins.pep.vs_ml2unfiltered.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c6.g634.i1	MLRB263543	33.7	323	191	11	10	320	37	348	1.7e-47	188.7
#Hcv1.1.c11.g551.i1	MLRB263543	43.6	326	166	5	47	356	34	357	3.3e-71	267.7
grep MLRB263549 Hcv1a1d20200325_model_proteins.pep.vs_ml2unfiltered.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c2.g345.i1	MLRB263549	81.8	303	55	0	71	373	339	641	3.4e-146	516.9
#Hcv1.1.c6.g622.i1	MLRB263549	68.8	304	94	1	64	367	339	641	6.3e-121	433.0
#Hcv1.1.c2.g700.i2	MLRB263549	32.1	252	168	2	69	319	335	584	2.1e-39	162.2
#Hcv1.1.c3.g629.i1	MLRB263549	33.8	201	128	2	111	310	9	205	2.1e-23	109.0
#Hcv1.1.c6.g172.i1	MLRB263549	30.6	170	107	3	130	299	1	159	5.2e-15	81.6
#Hcv1.1.c12.g332.i1	MLRB263549	26.2	237	132	2	126	319	18	254	9.0e-19	93.6
grep MLRB263551 Hcv1a1d20200325_model_proteins.pep.vs_ml2unfiltered.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c2.g345.i1	MLRB263551	82.8	379	65	0	25	403	21	399	6.5e-190	662.1
#Hcv1.1.c6.g622.i1	MLRB263551	70.2	393	116	1	5	397	8	399	3.5e-164	576.6
#Hcv1.1.c2.g700.i2	MLRB263551	29.0	366	258	2	33	397	26	390	3.2e-48	191.4
#Hcv1.1.c3.g629.i1	MLRB263551	30.3	399	269	6	1	394	1	395	8.6e-46	183.3
#Hcv1.1.c6.g172.i1	MLRB263551	27.4	332	232	4	96	426	68	391	6.1e-32	137.9
#Hcv1.1.c12.g332.i1	MLRB263551	29.6	378	251	6	31	399	26	397	3.2e-40	164.9






# for photoproteins
echo "# for photoproteins" >> Hcv1a1d20200325_biolum_results.txt
~/diamond-latest/diamond blastp -q hcal_nonpp.fa -d Hcv1a1d20200325_model_proteins.pep -o hcal_nonpp_from_txomes_vs_Hcv1.tab
cat hcal_nonpp_from_txomes_vs_Hcv1.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcal_nonPP	Hcv1.1.c9.g243.i1	76.5	196	43	1	1	193	1	196	2.2e-85	312.0
#Pbac_9913490_9953503_9853218	Hcv1.1.c9.g243.i1	88.1	185	22	0	9	193	12	196	1.7e-93	339.0
#Hormiphora_californensis_AOI27780.1	Hcv1.1.c13.g433.i1	66.8	241	73	3	1	239	1	236	2.0e-88	322.4
#Hormiphora_californensis_AOI27780.1	Hcv1.1.c13.g433.i2	66.8	241	73	3	1	239	1	236	3.4e-88	321.6
#Hcal_2OGFe1	Hcv1.1.c11.g551.i1	63.8	365	117	5	3	356	1	361	6.0e-129	457.6
#Hcal_2OGFe2	Hcv1.1.c6.g634.i1	75.8	335	80	1	12	345	5	339	1.4e-154	542.7
#Hcal_2OGFe2	Hcv1.1.c11.g551.i1	30.6	337	204	7	12	339	45	360	2.9e-40	162.9
#Dgla_FYY1	Hcv1.1.c11.g551.i1	42.0	345	175	7	16	356	33	356	3.4e-71	265.8
#Dgla_FYY1	Hcv1.1.c6.g634.i1	33.2	328	184	11	38	347	10	320	2.1e-49	193.4

echo "# for group 1" >> Hcv1a1d20200325_biolum_results.txt
# group 1, 6 inside this gene
grep ML085729a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c3.g1031.i2	ML085729a	55.6	288	123	2	137	421	22	307	1.7e-19	94.0

echo "# for group 2" >> Hcv1a1d20200325_biolum_results.txt
# group 2, sandwiched between these two
grep ML085714a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c3.g363.i1	ML085714a	33.2	398	230	10	19	392	381	766	1.6e-61	233.8
grep ML085716a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c3.g1150.i1	ML085716a	60.5	263	99	1	50	307	226	488	1.4e-80	296.6
#Hcv1.1.c3.g666.i1	ML085716a	37.0	257	153	3	73	324	226	478	6.5e-43	171.4
#Hcv1.1.c2.g1014.i1	ML085716a	33.7	273	139	8	83	332	239	492	6.4e-32	135.2
#Hcv1.1.c10.g206.i1	ML085716a	32.0	253	159	6	68	314	228	473	1.3e-27	120.9
#Hcv1.1.c11.g638.i2	ML085716a	30.0	297	175	8	48	339	213	481	6.4e-26	115.2
#Hcv1.1.c1.g467.i1	ML085716a	26.2	225	138	7	73	276	226	443	8.0e-08	54.7
#Hcv1.1.c2.g39.i1	ML085716a	29.0	290	180	7	76	364	230	494	5.1e-21	99.0
#Hcv1.1.c2.g864.i1	ML085716a	30.5	275	160	7	68	329	220	476	2.5e-19	93.2
#Hcv1.1.c2.g1255.i1	ML085716a	29.8	255	161	6	80	324	230	476	2.0e-18	90.1
#Hcv1.1.c5.g65.i1	ML085716a	27.4	266	151	8	66	319	230	465	2.1e-10	63.5
#Hcv1.1.c5.g309.i1	ML085716a	28.1	260	150	10	94	337	238	476	1.2e-13	74.3
#Hcv1.1.c5.g559.i1	ML085716a	23.1	346	157	10	63	398	230	476	6.9e-05	45.4
#Hcv1.1.c6.g46.i1	ML085716a	27.1	262	170	5	73	331	230	473	9.6e-17	84.7
#Hcv1.1.c7.g586.i1	ML085716a	28.5	165	103	4	59	220	226	378	3.1e-10	63.2
#Hcv1.1.c10.g263.i1	ML085716a	26.4	280	175	8	73	344	226	482	9.3e-12	68.2


echo "# for group 3" >> Hcv1a1d20200325_biolum_results.txt
# group 3, between these
grep ML215419a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c13.g148.i2	ML215419a	91.9	149	12	0	1	149	1	149	2.0e-71	267.3
grep ML215423a Hcv1a1d20200325_model_proteins.pep.vs_ml2.tab >> Hcv1a1d20200325_biolum_results.txt
#Hcv1.1.c10.g448.i2	ML215423a	64.2	430	147	4	99	522	1	429	1.5e-170	596.7
#Hcv1.1.c8.g472.i1	ML215423a	35.9	145	72	9	1269	1408	6	134	2.5e-10	65.5



################################################################################
################################################################################
################################################################################

# checking against photocyte proteins from Sebe-pedros 2018
# group 51

~/diamond-latest/diamond blastp -q sebepedros2018_photocyte_cluster_Mlei_prots.fasta -d Hcv1a1d20200325_model_proteins.pep.dmnd -o sebepedros2018_photocyte_cluster_vs_Hcv1a1d20200325.tab











