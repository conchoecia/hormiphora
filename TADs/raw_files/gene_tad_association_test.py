# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 13:02:13 2020

@author: Jakob McBroome
"""
#%%
import numpy as np
import seaborn as sns
import pandas as pd
import math
from scipy.stats import percentileofscore
#%%
predf = {k:[] for k in ['Chro','Start','Stop','Score']}
with open('../homer_tads/cteno_homer_1k_filtered_tads.bedgraph') as tadin:
    for entry in tadin:
        chro, start, stop, null, score = entry.strip().split()
        predf['Chro'].append(chro)
        predf['Start'].append(int(start))
        predf['Stop'].append(int(stop))
        predf['Score'].append(float(score))
datadf = pd.DataFrame(predf)
del predf
datadf
#%%
#snagged this from stackoverflow because lazy.
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="right")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]
#%%
lookup = {} #a simple dictionary structure will be easier.
with open('cteno_genes_pos.bed') as gin:
    for entry in gin:
        chro, start, stop = entry.strip().split()
        if chro not in lookup:
            lookup[chro] = []
        lookup[chro].append(int(start))
        lookup[chro].append(int(stop)) #in this case we don't particularly care about which goes with what.
for k,v in lookup.items():
    lookup[k] = np.array(v)
lookup
#%%
ns = []
ne = []
for i, d in datadf.iterrows():
    ns.append(abs(d.Start - find_nearest(lookup[d.Chro],d.Start)))
    ne.append(abs(d.Stop - find_nearest(lookup[d.Chro],d.Stop)))
datadf = datadf.assign(NearStart = ns)
datadf = datadf.assign(NearEnd = ne)
#%%
sns.scatterplot(x = (datadf.NearStart), y = (datadf.NearEnd))
#%%
#get the nearest chro-loc point set and save it, so it can be printed to a bed-like and added to the visualization
with open("nearest_gene_edges.txt", 'w+') as outf:
    for i, d in datadf.iterrows():
        sv = find_nearest(lookup[d.Chro],d.Start)
        ev = find_nearest(lookup[d.Chro],d.Stop)
        print(d.Chro + '\t' + str(sv) + '\t' + str(sv+10), file = outf)
        print(d.Chro + '\t' + str(ev) + '\t' + str(ev+10), file = outf)
#%%
sns.scatterplot(x = np.log10(datadf['NearStart']), y = np.log10(datadf['NearEnd']), data = datadf) 
#%%
sns.distplot(datadf['Stop']-datadf['Start'])
#%%
lengths = datadf['Stop']-datadf['Start']
length_nooutlier = [l for l in lengths if l < 150000]
print("Mean length without outliers: {:.4f}".format(np.mean(length_nooutlier)) )
print("Mean length, all data: {:.4f}".format(np.mean(lengths)))
print("Median length, no outliers: {:.4f}".format(np.median(length_nooutlier)))
print("Std err, no outliers: {:.4f}".format(np.std(length_nooutlier)))
#%%
clens = {}
with open('UCSC_Hcal_v1.fa.fai') as indexin:
    for entry in indexin:
        spent = entry.strip().split()
        chro = spent[0]
        clen = int(spent[1])
        clens[chro] = clen
print(clens)
#%%
def permuter(datadf, lookup, clens, pnum = 1000):
    #one thousand times, generate a random set of starts and ends with the same distribution, then measure the distance to the nearest gene edge
    #get that complete distribution of distances and rank the real one among them.
    distributions = {'starts':[], 'ends':[]}
    for p in range(pnum):
        pnstarts = []
        pnends = []
        chros = []
        starts = []
        ends = []
        for i, d in datadf[datadf['Stop'] - datadf['Start'] < 150000].iterrows():
            leng = d.Stop - d.Start
            rs = np.random.randint(0, clens[d.Chro] - leng - 1)
            re = rs + leng
            chros.append(d.Chro)
            starts.append(rs)
            ends.append(re)
        #proceed to finding nearest
        for i, c in enumerate(chros):
            pnstarts.append(abs(starts[i] - find_nearest(lookup[c],starts[i])))
            pnends.append(abs(ends[i] - find_nearest(lookup[c],ends[i])))
        distributions['starts'].append(pnstarts)
        distributions['ends'].append(pnends)
    return distributions
pdists = permuter(datadf, lookup, clens, pnum = 1000)
#%%
#convert these to mean values.
psmedian = [np.median(d) for d in pdists['starts']]
pemedian = [np.median(d) for d in pdists['ends']]
#%%
#sns.distplot(psmeans)
print(np.median(datadf[datadf['Stop'] - datadf['Start'] < 150000].NearStart))
print(np.percentile(psmedian, [5,25,50,75,95]))
print(percentileofscore(psmedian, np.median(datadf[datadf['Stop'] - datadf['Start'] < 150000].NearStart), kind='strict'))
print(percentileofscore(pemedian, np.median(datadf[datadf['Stop'] - datadf['Start'] < 150000].NearEnd), kind='strict'))
#pvalue = .2/100 = .002. 3.8/100 = .038.
print(percentileofscore(psmedian, np.median(datadf[datadf['Stop'] - datadf['Start'] < 150000][datadf.Score > 2].NearStart), kind='strict'))
print(percentileofscore(pemedian, np.median(datadf[datadf['Stop'] - datadf['Start'] < 150000][datadf.Score > 2].NearEnd), kind='strict'))
print(percentileofscore(psmedian + pemedian, np.median(datadf[datadf['Stop'] - datadf['Start'] < 150000].NearStart.tolist() + datadf[datadf['Stop'] - datadf['Start'] < 150000].NearEnd.tolist())))
#%%
#On request from Darrin, investigate the proportion of relative orientations for genes within TADs. 
#for the purpose of this, I'll define "within" as having any part of their span overlapping with the TAD.
lookup_stranded = {'pos':{},'neg':{}} #a simple dictionary structure will be easier.
with open('cteno_genes_pos.bed') as gin:
    for entry in gin:
        chro, start, stop = entry.strip().split()
        if chro not in lookup_stranded['pos']:
            lookup_stranded['pos'][chro] = []
        lookup_stranded['pos'][chro].append((int(start), int(stop)))
for k,v in lookup_stranded['pos'].items():
    lookup_stranded['pos'][k] = np.array(v)
with open('cteno_genes_neg.bed') as gin:
    for entry in gin:
        chro, start, stop = entry.strip().split()
        if chro not in lookup_stranded['neg']:
            lookup_stranded['neg'][chro] = []
        lookup_stranded['neg'][chro].append((int(start), int(stop)))
for k,v in lookup_stranded['neg'].items():
    lookup_stranded['neg'][k] = np.array(v)
#%%
pcount = 0
ncount = 0
for strand, scaff in lookup_stranded.items():
    for sc, d in scaff.items():
        if strand == 'pos':
            pcount += len(d)
        elif strand == 'neg':
            ncount += len(d)
print(pcount,ncount)
#%%
#now the intersection.
#two spans intersect if and only if both the start and the end of one are not both less than the start of the other or greater than the end of the other
#I feel like there must be an easier way to express that. Whatever. I'm not great with this math logic stuff a lot of the time. Sigh.
props = []
for i, d in datadf.iterrows():
    #stupid not smart algorithm because lazy and sleepy.
    pos = 0
    neg = 0
    for strand in ['pos','neg']:
        for gene in lookup_stranded[strand][d.Chro]:
            if gene[1] < d.Start:
                continue
            if gene[0] > d.Stop:
                break
            else:
                if strand == 'pos':
                    pos += 1
                elif strand == 'neg':
                    neg += 1
    props.append((pos, neg))
props
#%%
rp = [d[0]/sum(d) for d in props]
sns.distplot(rp)
#%%
from scipy.stats import fisher_exact
for p in props:
    r,pv = fisher_exact([p,[pcount,ncount]])
    if pv < .05/len(props):
        print(p,pv)
#%%
sns.distplot([sum(d) for d in props])












#%%
#there's somewhat of a signal of existing near it. What about interrupting it? Can go back to the start/end and check that.
#god I feel like I've been here before.
lookup2 = {} #a simple dictionary structure will be easier.
with open('cteno_genes.bed') as gin:
    for entry in gin:
        chro, start, stop = entry.strip().split()
        if chro not in lookup2:
            lookup2[chro] = []
        lookup2[chro].append((int(start),int(stop)))
lookup2
#%%
startbreaks_real = []
endbreaks_real = []
for i, d in datadf.iterrows():
    #using a stupid algorithm because I'm lazy and don't particularly need it to run fast
    covers_start = 0
    covers_end = 0
    for gene in lookup2[d.Chro]:
        if gene[1] < d.Start:
            continue #keep looking
        elif gene[1] >= d.Start and gene[0] < d.Start:
            covers_start += 1 #if the gene starts before it starts and ends after
        elif gene[0] <= d.Stop and gene[1] > d.Stop:
            covers_end += 1 #if the gene starts before it ends and ends after
        elif gene[0] > d.Stop:
            break #go to the next TAD.
    startbreaks_real.append(covers_start)
    endbreaks_real.append(covers_end)
print(startbreaks_real,endbreaks_real)
#%%
def permuter2(datadf, lookup, clens, pnum = 1000):
    #one thousand times, generate a random set of starts and ends with the same distribution, then measure the distance to the nearest gene edge
    #get that complete distribution of distances and rank the real one among them.
    breakage_sets = []
    for p in range(pnum):
        startbreaks = []
        endbreaks = []
        chros = []
        starts = []
        ends = []
        for i, d in datadf.iterrows():
            leng = d.Stop - d.Start
            rs = np.random.randint(0, clens[d.Chro] - leng - 1)
            re = rs + leng
            chros.append(d.Chro)
            starts.append(rs)
            ends.append(re)
        #proceed to finding nearest
        for i, c in enumerate(chros):
            #using a stupid algorithm because I'm lazy and don't particularly need it to run fast
            covers_start = 0
            covers_end = 0
            for gene in lookup2[c]:
                if gene[1] < starts[i]:
                    continue #keep looking
                elif gene[1] >= starts[i] and gene[0] < starts[i]:
                    covers_start += 1 #if the gene starts before it starts and ends after
                elif gene[0] <= ends[i] and gene[1] > ends[i]:
                    covers_end += 1 #if the gene starts before it ends and ends after
                elif gene[0] > ends[i]:
                    break #go to the next TAD.
            startbreaks.append(covers_start)
            endbreaks.append(covers_end)
            
        breakage_sets.append((startbreaks,endbreaks))
    return breakage_sets
pdists = permuter2(datadf, lookup, clens, pnum = 1000)
#%%
#get the mean number of tad boundaries which intersect a gene span across pdists.
intersects = [(len([v for v in p[0] if v != 0]), len([v for v in p[1] if v != 0])) for p in pdists]
print(np.mean(intersects))
rsn = len([v for v in startbreaks_real if v != 0])
esn = len([v for v in endbreaks_real if v != 0])
allpts = np.array(intersects).flatten()
print(rsn, esn)
print(np.percentile(allpts, [5,25,50,75,95]))
