__author__ = 'ranko'

from RBioTools import loadAllFastas
from RBioTools import hashFastas
from RBioTools import getFaSeq
from RBioTools import getFaHeader
from RBioTools import saveAllFastas
from RBioTools import countAllAAs
from RBioTools import countDiMers
from RBioTools import countEGroups
from RBioTools import countDiEGroup
from RBioTools import munchSpaces

import os

# ------------------------------------
# misc: parses BLAST tab output into 'pretty' html table
# must be in format:
# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, % subject coverage, % positives, query length, subject length
#
#
# ------------------------------------


# ----------------------------------
# scored ToxBits vectorizer
#  for generating positive training set,classifyNoTB should be set to 0,
#  otherwise it should be set to 1
#
# ----------------------------------
def scoredToxBitsVectorize(inp,toxBits,hmmer, classifyNoTB=1,verbose=False, tmpfolder = '',repor=-1,trimHdr=0):

    if tmpfolder == '':
        tmpfolder = './'
    elif not tmpfolder[-1] == '\/':
        tmpfolder = tmpfolder+'/'
# ---- STEP 0: PREPROCESS TOXBITS ----
    if verbose:
        print ' scored Toxbits vectorizing started ... '
        print ' -------------------------------------------------------------------------- '
        print '    0) analyzing toxBits'
    allBits = []
    with open(toxBits) as iF:
        for l in iF:
            if l[0:4] == 'NAME':
                nm = munchSpaces(l).split(' ')[1].strip()
                allBits.append(nm)

    allBits = sorted(allBits)
    if verbose:
        print '    --> DONE; Found ',len(allBits),'tox-bits'

    # ---- GO THROUHG FASTAS 1 by 1:
    allResults = {} # hdr -> [ (toxBit1,score1), ... (toxBitn, scoren) ]
    if verbose:
        print '    -> HMMERING ...'
    saveAllFastas(inp,tmpfolder+'__tmp_tbv_input.fa')
    hmms = hmmer+'hmmsearch -E 1.0e+2 --cpu 8 --notextw --noali --tblout '+tmpfolder+'__tmp_tbv_res.hmmtab '+toxBits+' '+tmpfolder+'__tmp_tbv_input.fa > '+tmpfolder+'__tmp_tbv_res.out'
    if verbose:
        print hmms
    os.system(hmms)
    os.system('wait')
    # --- prep all headers (not all with have hits!)
    allHdrs = set()
    allFa = loadAllFastas(tmpfolder+'__tmp_tbv_input.fa')
    for f in allFa:
        allHdrs.add(getFaHeader(f).strip())

    # --- STEP 2: PARSE
    if verbose:
        print '    -> PARSING ... '
    fndToxBits = {} # -> hdr -> toxbits
    with open(tmpfolder+'__tmp_tbv_res.hmmtab') as iF:
        for l in iF:
            if not l[0] == '#':
                ls = munchSpaces(l.strip()).split(' ')
                #print ls
                if ls[18].strip().replace(' ','') == '-':
                    hdr = '>'+ls[0]
                else:
                    hdr = '>'+ls[0]+' '+' '.join(ls[18:])
                fndToxBits = ( (ls[2],ls[5]) )
                try:
                    allResults[hdr].append(fndToxBits)
                except:
                    allResults[hdr] = [fndToxBits]

    # -- check consistency:
    if verbose:
        print '        -> CONSISTENCY CHECK: nothing should be displayed!'
        good = True
        #d print allHdrs
        for h in allResults.keys():
            #d print h
            if h not in allHdrs:
                good = False
                print ' WARNING : ',h,'not in headers, bug is lurking somewhere...'
        if good:
            print '         -> CONSISTENCY GOOD (Gog like) !'
        print '        -> DONE PARSING ... '

    # ---- VECTORIZING -----
    if verbose:
        print '    -> VECTORIZING ...'
    allVecs = []
    title = ['qID']
    allVecs.append(title)
    for t in allBits:
        title.append(t)
    fn = 0
    p = 0
    rm = len(allHdrs)
    cnothing = 0
    hdrsWTB = set()
    for k in allHdrs:
        nothing = True
        fn+=1
        p+=1
        if p >= rm/100*5:
            if verbose:
                print '   vectorizing',fn,'/',rm,'[',int(float(fn)/float(rm)*100.0),'% ]'
        p = 0
        #print k
        oneVec = []
        if trimHdr == 1:
            if ' ' in k:
                oneVec.append(k.split(' ')[0])
            else:
                oneVec.append(k)
        else:
            oneVec.append(k)
        try:
            # sort toxbits
            sortb = sorted(allResults[k],key = lambda k : k[0])
            #print sortb
            #if saveWTB == '1':
            #    hdrsWTB.add(k)
        except:
            #if args.classifyNoTB == '1':
            sortb = [(0,0)]
            #print ' warning: no toxbits for ',k
    #print sortb
        ctb = 0
        for t in allBits:
            fnd = False
            for tk in sortb:
                if str(t) == str(tk[0]):
                    oneVec.append(tk[1])
                    ctb+=1
                    #print t,tk[0],tk[1]
                    fnd = True
                    nothing = False
                    break
            if not fnd:
                oneVec.append('0.0')
                ctb+=1
        #if (args.classifyNoTB == '1' and nothing) or fnd:
        if (nothing and not classifyNoTB == '1'):
            cnothing +=1
        else:
            allVecs.append(oneVec)
        if verbose:
            if (len(oneVec) != len(allBits)+1) or ctb != len(allBits):
                print 'Warning:', k,'bits:',ctb,'vec lt:',len(oneVec),'(should be',ctb+1,')'
    if verbose:
        print 'Scored ToxBits vectorizer done:'
        print ' -> vectorized total of ',fn,'seqs'
        print ' -> sequences with 0 toxbits: ',cnothing
        print ' -> input had ',len(allHdrs),'seqs'
        print ' -> saved ',len(allVecs),'vectors'

    return allVecs

# ----------------------------------
# biFreq vectorizer
# ----------------------------------
def biFreqVectorize(inp,verbose=False, inputIsFile = False, trimHeader=0):
    tit = ['qID','qLT']
    hesh = {}
    aas = 'ACDEFGHIKLMNPQRSTVWYX'
    for a in aas:
        tit.append(a)
    for a1 in aas:
        for a2 in aas:
            hesh[str(str(a1)+str(a2))] = 0.0
    for h in sorted(hesh.keys()):
        tit.append(h)
    if verbose:
        print ' -> bifreq vectorizing... '
        print ' -------------------------------------------------------------------------- '
    allvecs = []
    allvecs.append(tit)

    if inputIsFile:
        iF = loadAllFastas(inp)
    else:
        iF = inp

    for f in iF:
        hdr = getFaHeader(f).strip()
        if trimHeader == 1:
            hdr = getFaHeader(f).strip().split(' ')[0].strip()

        onevec = [hdr,len(getFaSeq(f).strip())]
        onemers = countAllAAs(f)
        for a in sorted(onemers.keys()):
            onevec.append(onemers[a])
            #print len(onemers.keys())
        dimers = countDiMers(f)
        for a in sorted(dimers.keys()):
            #print len(dimers.keys())
            onevec.append(dimers[a])

        if len(onevec) != 464:
            print ' -> BIFREQ ERROR:',getFaHeader(f),'vector l = ',len(onevec),'should be 464!'
            #print onevec
            exit()
        allvecs.append(onevec)

    if verbose:
        print ' -> bifreq vectorizing done, processed ',len(allvecs)-1,'sequences!'
        print ' -------------------------------------------------------------------------- '

    return allvecs

# ----------------------------------
# threeBLAST vectorizer
# ----------------------------------

def triBlastEnhVectorize (iSeq,dbNeg1='',dbNeg2='',dbPos='',cutoff=1.0e+1,blastPath='',verbose=False, inputIsFile = False, tmpfolder='', trimHeader = 0):
    if tmpfolder == '':
        tmpfolder = './'
    elif not tmpfolder[-1] == '\/':
        tmpfolder = tmpfolder+'/'
    #print iSeq
    c = 0
    res = []
    if verbose:
        print ' -> triBlastEnhanced vectorizing... '
        print ' -------------------------------------------------------------------------- '
    vecT = ['qID','qLT','Tox_sLT','Tox_qLT/sLT','Tox_bScore','Tox_qCov','Tox_pIdent','Tox_pPos','C_SP_sLT','C_SP_qLT/sLT','C_SP_bScore','C_SP_qCov','C_SP_pIdent','C_SP_pPos','C_TR_sLT','C_TR_qLT/sLT','C_TR_bScore','C_TR_qCov','C_TR_pIdent','C_TR_pPos', \
            'percA','percC','percD','percE','percF','percG','percH','percI','percK','percL','percM','percN','percP','percQ','percR','percS','percT','percV','percW','percY','percX']
    if verbose:
        print '    -> BLASTING vs negative db 1 '
    # blast them
    saveAllFastas(iSeq,tmpfolder+'__tmpfile.fa')
    # blast 1 (neg1)
    bCmd = blastPath+' -max_target_seqs 10 -query '+tmpfolder+'__tmpfile.fa -outfmt "7 std qcovs ppos qlen slen stitle" -evalue '+str(cutoff)+' -db '+dbNeg1+' -num_threads 20 > '+tmpfolder+'__tb__tmp_neg1.bres'
    if verbose:
        print bCmd
    os.system(bCmd)
    # blast 2 (neg2)
    if verbose:
        print '    -> BLASTING vs negative db 2 '
    bCmd = blastPath+' -max_target_seqs 10 -query '+tmpfolder+'__tmpfile.fa -outfmt "7 std qcovs ppos qlen slen stitle" -evalue '+str(cutoff)+' -db '+dbNeg2+' -num_threads 20 > '+tmpfolder+'__tb__tmp_neg2.bres'
    if verbose:
        print bCmd
    os.system(bCmd)
    # blast 3 (pos)
    if verbose:
        print '    -> BLASTING vs positive DB '
    bCmd = blastPath+' -max_target_seqs 10 -query '+tmpfolder+'__tmpfile.fa -outfmt "7 std qcovs ppos qlen slen stitle" -evalue '+str(cutoff)+' -db '+dbPos+' -num_threads 20 > '+tmpfolder+'__tb__tmp_pos.bres'
    if verbose:
        print bCmd
    os.system(bCmd)
    # now load all three
    allHdrs = set()
    afh = loadAllFastas(tmpfolder+'__tmpfile.fa')
    faHesh = hashFastas(afh)
    for f in afh:
        allHdrs.add(getFaHeader(f).strip())
    #print allHdrs
    # PARSE IT !!!
    if verbose:
        print '    -> vectorizing BLAST results '
    # hashes: header -> statistic

    toxSL = {}
    toxQLSL = {}
    toxBS = {}
    toxQC = {}
    toxPI = {}
    toxPP = {}
    cspSL = {}
    cspQLSL = {}
    cspBS = {}
    cspQC = {}
    cspPI = {}
    cspPP = {}
    ctrSL = {}
    ctrQLSL = {}
    ctrBS = {}
    ctrQC = {}
    ctrPI = {}
    ctrPP = {}
    # construct vector
    # -> basically: after all hashes are full, go through all headers
    # and construct vector for appropriate header by grabbing from hashes
    # there will be empty stuff, in which case just put 0 in there
    with open(tmpfolder+'__tb__tmp_pos.bres') as iF:
        for l in iF:
            l = l.strip()
            if not l[0] == '#':
                if l.split('\t')[0] not in toxSL.keys():
                    # print l.split('\t')[0]
                    toxSL[l.split('\t')[0]] = float(l.split('\t')[15])
                    toxBS[l.split('\t')[0]] = float(l.split('\t')[11])
                    toxQC[l.split('\t')[0]] = float(l.split('\t')[12])/100.0
                    toxPI[l.split('\t')[0]] = float(l.split('\t')[2])/100.0
                    toxPP[l.split('\t')[0]] = float(l.split('\t')[13])/100.0
    with open(tmpfolder+'__tb__tmp_neg1.bres') as iF:
        for l in iF:
            l = l.strip()
            if not l[0] == '#':
                if l.split('\t')[0] not in cspSL.keys():
                    cspSL[l.split('\t')[0]] = float(l.split('\t')[15])
                    cspBS[l.split('\t')[0]] = float(l.split('\t')[11])
                    cspQC[l.split('\t')[0]] = float(l.split('\t')[12])/100.0
                    cspPI[l.split('\t')[0]] = float(l.split('\t')[2])/100.0
                    cspPP[l.split('\t')[0]] = float(l.split('\t')[13])/100.0
    with open(tmpfolder+'__tb__tmp_neg2.bres') as iF:
        for l in iF:
            l = l.strip()
            if not l[0] == '#':
                if l.split('\t')[0] not in ctrSL.keys():
                    ctrSL[l.split('\t')[0]] = float(l.split('\t')[15])
                    ctrBS[l.split('\t')[0]] = float(l.split('\t')[11])
                    ctrQC[l.split('\t')[0]] = float(l.split('\t')[12])/100.0
                    ctrPI[l.split('\t')[0]] = float(l.split('\t')[2])/100.0
                    ctrPP[l.split('\t')[0]] = float(l.split('\t')[13])/100.0
    #vecT = ['qID','qLT','Tox_sLT','Tox_qLT/sLT','Tox_bScore','Tox_qCov','Tox_pIdent','Tox_pPos',\
    #        'C_SP_sLT','C_SP_qLT/sLT','C_SP_bScore','C_SP_qCov','C_SP_pIdent','C_SP_pPos',\
    #        'C_TR_sLT','C_TR_qLT/sLT','C_TR_bScore','C_TR_qCov','C_TR_pIdent','C_TR_pPos']

    # make all vectors
    e1added = False
    if verbose:
        print '    -> calculating Amino Acid frequencies... '

    for hdr in allHdrs:
        htr = hdr.split(' ')[0][1:]
        if htr not in toxSL.keys() and htr not in cspSL.keys() and htr not in ctrSL.keys():
            print ' warning: ',htr,' has no close seq in any test DB'
        #print htr
        seq = faHesh[hdr.strip()].strip()
        vec = []
        sQL = len(seq)
        #  -------- toxin query LT / subject LT
        try:
            sToxQLSL = float(sQL)/max(1,toxSL[htr])
        except:
            #print ' warning',htr,'not in toxSL'
            sToxQLSL = 0.0
        #  -------- neg1 query LT / subject LT
        try:
            sCspQLSL = float(sQL)/max(1,cspSL[htr])
        except:
            #print ' warning',htr,'not in cspSL'
            sCspQLSL = 0.0
        #  -------- neg2 query LT / subject LT
        try:
            sCtrQLSL = float(sQL)/max(1,ctrSL[htr])
        except:
            #print ' warning',htr,'not in ctrSL'
            sCtrQLSL = 0.0
        # -------------- toxin subject lt ----------
        try:
            sToxSL = toxSL[htr]
        except:
            #print ' warning: ',htr,'not in toxSL'
            sToxSL = 0.0
        # -------- sToxBS (toxin bit score) -------
        try:
            sToxBS = toxBS[htr]
        except:
            #print ' warning: ',htr,'not in toxBS'
            sToxBS = 0.0
        # -------- sToxQC (toxin query coverage) ------
        try:
            sToxQC = toxQC[htr]
        except:
            #print ' warning: ',htr,'not in sToxQC'
            sToxQC = 0.0
        # --------- sToxPI (toxin identity)
        try:
            sToxPI = toxPI[htr]
        except:
            #print ' warning: ',htr,'not in sToxPI'
            sToxPI = 0.0
        # --------- sToxPP (toxin percentage partial match) --------
        try:
            sToxPP = toxPP[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sToxPP = 0.0
        # --------- sCspSL (neg 1 subject lt) --------
        try:
            sCspSL = cspSL[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCspSL = 0.0
        # --------- sCspBS (neg 1 subject score) --------
        try:
            sCspBS = cspBS[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCspBS = 0.0
        # --------- sCspBS (neg 1 subject coverage) --------
        try:
            sCspQC = cspQC[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCspQC = 0.0
        # --------- sCspPI (neg 1) --------
        try:
            sCspPI = cspPI[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCspPI = 0.0
        # --------- sCspPP (neg 1) --------
        try:
            sCspPP = cspPP[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCspPP = 0.0
        # --------- sCspPP (neg 1) --------
        try:
            sCtrSL = ctrSL[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCtrSL = 0.0
        # --------- sCtrBS (neg 2) --------
        try:
            sCtrBS = ctrBS[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCtrBS = 0.0
        # --------- sCtrQC (neg 2) --------
        try:
            sCtrQC = ctrQC[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCtrQC = 0.0
        # --------- sCtrPI (neg 2) --------
        try:
            sCtrPI = ctrPI[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCtrPI = 0.0
        # --------- sCtrPI (neg 2) --------
        try:
            sCtrPP = ctrPP[htr]
        except:
            #print ' warning: ',htr,'not in sToxPP'
            sCtrPP = 0.0

        # ------- COUNT FEQS (1)
        vec = ['>'+htr,sQL,sToxSL,sToxQLSL,sToxBS,sToxQC,sToxPI,sToxPP,sCspSL,sCspQLSL,sCspBS,sCspQC,sCspPI,sCspPP,sCtrSL,sCtrQLSL,sCtrBS,sCtrQC,sCtrPI,sCtrPP]
        aaCounts = []
        aaCountsDick = countAllAAs(seq)
        for k in sorted(aaCountsDick.keys()):
            aaCounts.append(aaCountsDick[k])
        vec = vec + aaCounts
        res.append(vec)
        # ------- COUNT E groups (1)
        r = countEGroups(seq)
        #print r
        rk = sorted(r.keys())
        #print rk
        for k in rk:
            if not e1added:
                vecT.append(k)
            vec.append(r[k])
        # ------ COUNT AAFreq (2-mer)
        dimers = countDiMers(seq)
        for a in sorted(dimers.keys()):
            #print len(dimers.keys())
            vec.append(dimers[a])
            if not e1added:
                vecT.append(a)
        # ------ COUNT E groups (2-mer)
        dimers = countDiEGroup(seq)
        for a in sorted(dimers.keys()):
            #print len(dimers.keys())
            vec.append(dimers[a])
            if not e1added:
                vecT.append(a)
        e1added = True
        if len(res) % 1000 == 0 and verbose:
            print '    -> ',len(res),'seqs done!'
    #print sorted(aaCountsDick.keys())
    #print aaCounts
    if verbose:
        print ' -> triBlastEnhanced vectorization done, processed ',len(res),'sequences'
        print ' -------------------------------------------------------------------------- '

    res = [vecT]+res
    return res
# -----------------------------------------------------------------
