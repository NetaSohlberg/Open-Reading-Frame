'''
Open Reading Frame (ORF)/ Neta Sohlberg

"In molecular genetics, an open reading frame (ORF) is the part of a reading frame that has the ability to be translated.
 An ORF is a continuous stretch of codons that may begin with a start codon (usually ATG) and ends at a stop codon (usually TAA, TAG or TGA)" (Wikipedia)
 
The program finds several possible ORF in DNA sequence and its translation into protein
'''

#the function return a frame according to the num it get (+1,+2,+3,-1,-2 or -3)
def find_frame(seq,num):
    newseq=""
    #if num<0- reverse the seq
    if num<0:
        newseq=reverse(seq)
        num=-num
    else:
        nuc=['A','C','T','G']
        for n in seq:
            if n in nuc:
                newseq+=n
    #frame is a list of threesome
    #count=[i for i in range(len(newseq)-3)]
    frame=[newseq[0:num-1]]
    i=num-1
    while i<(len(newseq)-3):
        frame+=[newseq[i:i+3]]
        i+=3
    #call the function "find_ORF" to print the ORF in the frame
    find_ORF(frame)
#end function

#reverse a seq
def reverse(seq):
    newseq=""
    dic={'A':'T','C':'G','T':'A','G':'C'}
    for nuc in seq:
        if nuc in dic:
            newseq+=dic[nuc]
    return newseq[::-1]
#end function

#find and print ORF and its translation, using the function "print_proteins"
def find_ORF(frame):
    if (len(frame))<100:
        return
    start=['ATG']
    stop=['TAG', 'TGA', 'TAA']
    #chek if there is start codon and stop codon
    if bool(set(start).intersection(frame))and bool(set(stop).intersection(frame)):
        i1=frame.index('ATG')
        for i in frame:
            if i in stop:
                i2=frame.index(i)
                break
        #check if the length is 300+, if not- try to find another ORF
        if (i2*3-i1*3)<100:
            find_ORF(frame[i2+1:len(frame)])
        else:
            #print details
            print ("start:",i1*3+1,end=' ')
            print ("stop:",i2*3+1,end=' ')
            print ("length:",i2*3-i1*3+1,end=' ')
            #printfirst and last nucleotides
            print (frame[i1],frame[i1+1],"\t",frame[i2-1],frame[i2])
            #print ORF sequence
            print ("ORF:",frame[i1:i2+1])
            #print transletion
            print_proteins(frame[i1:i2+1])
    else:
        print ("no ORF")
#end function

#print the translation
def print_proteins(ORF_frame):
    code = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"s", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
    "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G"}
    print ("translation:")
    for codon in ORF_frame:
        print (code[codon],end=' ')
#end function

#main
f=open("seq.txt","r")
f.readline()
seq=""
for line in f:
    seq+=line
#seq=f.read()

for i in [-3,-2,-1,1,2,3]:
    print ("\n\nFRAME",i)
    find_frame(seq.upper(),i)

f.close()
#end main


'''
#output for Cloning vector pACYC184
FRAME -3
start: 82 stop: 430 length: 349 ATG ACC 	ATA TAG
ORF: ['ATG', 'ACC', 'CAG', 'AGC', 'GCT', 'GCC', 'GGC', 'ACC', 'TGT', 'CCT', 'ACG', 'AGT', 'TGC', 'ATG', 'ATA', 'AAG', 'AAG', 'ACA', 'GTC', 'ATA', 'AGT', 'GCG', 'GCG', 'ACG', 'ATA', 'GTC', 'ATG', 'CCC', 'CGC', 'GCC', 'CAC', 'CGG', 'AAG', 'GAG', 'CTG', 'ACT', 'GGG', 'TTG', 'AAG', 'GCT', 'CTC', 'AAG', 'GGC', 'ATC', 'GGT', 'CGA', 'CGC', 'TCT', 'CCC', 'TTA', 'TGC', 'GAC', 'TCC', 'TGC', 'ATT', 'AGG', 'AAG', 'CAG', 'CCC', 'AGT', 'AGT', 'AGG', 'TTG', 'AGG', 'CCG', 'TTG', 'AGC', 'ACC', 'GCC', 'GCC', 'GCA', 'AGG', 'AAT', 'GGT', 'GCA', 'TGC', 'AAG', 'GAG', 'ATG', 'GCG', 'CCC', 'AAC', 'AGT', 'CCC', 'CCG', 'GCC', 'ACG', 'GGG', 'CCT', 'GCC', 'ACC', 'ATA', 'CCC', 'ACG', 'CCG', 'AAA', 'CAA', 'GCG', 'CTC', 'ATG', 'AGC', 'CCG', 'AAG', 'TGG', 'CGA', 'GCC', 'CGA', 'TCT', 'TCC', 'CCA', 'TCG', 'GTG', 'ATG', 'TCG', 'GCG', 'ATA', 'TAG']
translation:
M T Q S A A G T C P T S C M I K K T V I S A A T I V M P R A H R K E L T G L K A L K G I G R R S P L C D s C I R K Q P S S R L R P L S T A A A R N G A C K E M A P N S P P A T G P A T I P T P K Q A L M S P K W R A R S s P S V M S A I STOP 

FRAME -2
start: 37 stop: 202 length: 166 ATG AGC 	GGG TGA
ORF: ['ATG', 'AGC', 'AAA', 'CTG', 'AAA', 'CGT', 'TTT', 'CAT', 'CGC', 'TCT', 'GGA', 'GTG', 'AAT', 'ACC', 'ACG', 'ACG', 'ATT', 'TCC', 'GGC', 'AGT', 'TTC', 'TAC', 'ACA', 'TAT', 'ATT', 'CGC', 'AAG', 'ATG', 'TGG', 'CGT', 'GTT', 'ACG', 'GTG', 'AAA', 'ACC', 'TGG', 'CCT', 'ATT', 'TCC', 'CTA', 'AAG', 'GGT', 'TTA', 'TTG', 'AGA', 'ATA', 'TGT', 'TTT', 'TCG', 'TCT', 'CAG', 'CCA', 'ATC', 'CCT', 'GGG', 'TGA']
translation:
M S K L K R F H R S G V N T T T I s G S F Y T Y I R K M W R V T V K T W P I s L K G L L R I C F S S Q P I P G STOP 

FRAME -1
start: 7 stop: 442 length: 436 ATG GCA 	GCG TAA
ORF: ['ATG', 'GCA', 'ATG', 'AAA', 'GAC', 'GGT', 'GAG', 'CTG', 'GTG', 'ATA', 'TGG', 'GAT', 'AGT', 'GTT', 'CAC', 'CCT', 'TGT', 'TAC', 'ACC', 'GTT', 'TTC', 'CAT', 'GAG', 'CAA', 'ACT', 'GAA', 'ACG', 'TTT', 'TCA', 'TCG', 'CTC', 'TGG', 'AGT', 'GAA', 'TAC', 'CAC', 'GAC', 'GAT', 'TTC', 'CGG', 'CAG', 'TTT', 'CTA', 'CAC', 'ATA', 'TAT', 'TCG', 'CAA', 'GAT', 'GTG', 'GCG', 'TGT', 'TAC', 'GGT', 'GAA', 'AAC', 'CTG', 'GCC', 'TAT', 'TTC', 'CCT', 'AAA', 'GGG', 'TTT', 'ATT', 'GAG', 'AAT', 'ATG', 'TTT', 'TTC', 'GTC', 'TCA', 'GCC', 'AAT', 'CCC', 'TGG', 'GTG', 'AGT', 'TTC', 'ACC', 'AGT', 'TTT', 'GAT', 'TTA', 'AAC', 'GTG', 'GCC', 'AAT', 'ATG', 'GAC', 'AAC', 'TTC', 'TTC', 'GCC', 'CCC', 'GTT', 'TTC', 'ACC', 'ATG', 'GGC', 'AAA', 'TAT', 'TAT', 'ACG', 'CAA', 'GGC', 'GAC', 'AAG', 'GTG', 'CTG', 'ATG', 'CCG', 'CTG', 'GCG', 'ATT', 'CAG', 'GTT', 'CAT', 'CAT', 'GCC', 'GTC', 'TGT', 'GAT', 'GGC', 'TTC', 'CAT', 'GTC', 'GGC', 'AGA', 'ATG', 'CTT', 'AAT', 'GAA', 'TTA', 'CAA', 'CAG', 'TAC', 'TGC', 'GAT', 'GAG', 'TGG', 'CAG', 'GGC', 'GGG', 'GCG', 'TAA']
translation:
M A M K D G E L V I W D S V H P C Y T V F H E Q T E T F S S L W S E Y H D D F R Q F L H I Y S Q D V A C Y G E N L A Y F P K G F I E N M F F V S A N P W V S F T S F D L N V A N M D N F F A P V F T M G K Y Y T Q G D K V L M P L A I Q V H H A V C D G F H V G R M L N E L Q Q Y C D E W Q G G A STOP 

FRAME 1
start: 13 stop: 169 length: 157 ATG GTG 	GGG TGA
ORF: ['ATG', 'GTG', 'AAA', 'GTT', 'GGA', 'ACC', 'TCT', 'TAC', 'GTG', 'CCG', 'ATC', 'AAC', 'GTC', 'TCA', 'TTT', 'TCG', 'CCA', 'AAA', 'GTT', 'GGC', 'CCA', 'GGG', 'CTT', 'CCC', 'GGT', 'ATC', 'AAC', 'AGG', 'GAC', 'ACC', 'AGG', 'ATT', 'TAT', 'TTA', 'TTC', 'TGC', 'GAA', 'GTG', 'ATC', 'TTC', 'CGT', 'CAC', 'AGG', 'TAT', 'TTA', 'TTC', 'GGC', 'GCA', 'AAG', 'TGC', 'GTC', 'GGG', 'TGA']
translation:
M V K V G T S Y V P I N V S F S P K V G P G L P G I N R D T R I Y L F C E V I F R H R Y L F G A K C V G STOP 

FRAME 2
start: 163 stop: 544 length: 382 ATG CGC 	GCA TAA
ORF: ['ATG', 'CGC', 'ACC', 'CGT', 'TCT', 'CGG', 'AGC', 'ACT', 'GTC', 'CGA', 'CCG', 'CTT', 'TGG', 'CCG', 'CCG', 'CCC', 'AGT', 'CCT', 'GCT', 'CGC', 'TTC', 'GCT', 'ACT', 'TGG', 'AGC', 'CAC', 'TAT', 'CGA', 'CTA', 'CGC', 'GAT', 'CAT', 'GGC', 'GAC', 'CAC', 'ACC', 'CGT', 'CCT', 'GTG', 'GAT', 'CCT', 'CTA', 'CGC', 'CGG', 'ACG', 'CAT', 'CGT', 'GGC', 'CGG', 'CAT', 'CAC', 'CGG', 'CGC', 'CAC', 'AGG', 'TGC', 'GGT', 'TGC', 'TGG', 'CGC', 'CTA', 'TAT', 'CGC', 'CGA', 'CAT', 'CAC', 'CGA', 'TGG', 'GGA', 'AGA', 'TCG', 'GGC', 'TCG', 'CCA', 'CTT', 'CGG', 'GCT', 'CAT', 'GAG', 'CGC', 'TTG', 'TTT', 'CGG', 'CGT', 'GGG', 'TAT', 'GGT', 'GGC', 'AGG', 'CCC', 'CGT', 'GGC', 'CGG', 'GGG', 'ACT', 'GTT', 'GGG', 'CGC', 'CAT', 'CTC', 'CTT', 'GCA', 'TGC', 'ACC', 'ATT', 'CCT', 'TGC', 'GGC', 'GGC', 'GGT', 'GCT', 'CAA', 'CGG', 'CCT', 'CAA', 'CCT', 'ACT', 'ACT', 'GGG', 'CTG', 'CTT', 'CCT', 'AAT', 'GCA', 'GGA', 'GTC', 'GCA', 'TAA']
translation:
M R T R S R S T V R P L W P P P S P A R F A T W S H Y R L R D H G D H T R P V D P L R R T H R G R H H R R H R C G C W R L Y R R H H R W G R S G S P L R A H E R L F R R G Y G G R P R G R G T V G R H L L A C T I P C G G G A Q R P Q P T T G L L P N A G V A STOP 

FRAME 3
start: 22 stop: 136 length: 115 ATG TTG 	CAC TGA
ORF: ['ATG', 'TTG', 'GCA', 'CTG', 'ATG', 'AGG', 'GTG', 'TCA', 'GTG', 'AAG', 'TGC', 'TTC', 'ATG', 'TGG', 'CAG', 'GAG', 'AAA', 'AAA', 'GGC', 'TGC', 'ACC', 'GGT', 'GCG', 'TCA', 'GCA', 'GAA', 'TAT', 'GTG', 'ATA', 'CAG', 'GAT', 'ATA', 'TTC', 'CGC', 'TTC', 'CTC', 'GCT', 'CAC', 'TGA']
translation:
M L A L M R V S V K C F M W Q E K K G C T G A S A E Y V I Q D I F R F L A H STOP
'''

