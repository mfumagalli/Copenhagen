
def convertSyms(read,site):
    ''' read must be a Reads type object and site a Sites type object'''
    bases=""
    indexDelN=[]

    i=0
    while i <= len(read.base)-1:
        if read.base[i] in ['.',',']: #reference
            bases=str(bases+site.reference)
        elif read.base[i] in ['A','C','G','T','a','c','g','t']: #alternate
            bases=str(bases+read.base[i].upper())
        elif read.base[i] in ['^']: # start read, index +1 since the following character is mapping quality
            i=i+1
        elif read.base[i] in ['$']: #end read
            i=i #do nothing but keep for clarity
        elif read.base[i] in ['*','N','n','>','<']:  # asterisk is deleted base, N or n is undefined, > and < are reference skips, these will be then filtered out later on
            bases=str(bases+'X')
            indexDelN.append(i)
        elif read.base[i] in ['-','+']: # indel, skip to the the next non-indel base
            lenIndel = int(read.base[i+1])
            try: # Indels longer than 9
                lenIndel=int(str(lenIndel)+str(int(read.base[i+2])))+1
            except:
                pass
            i=i+lenIndel+1
        i=i+1
    return(bases,indexDelN)

def calcNonMajorCounts(read):

    alleles =['A','C','G','T']
    counts = [0,0,0,0]

    if len(read.base)>0:

        for i in range(len(read.base)):
            counts[alleles.index(read.base[i])]+=1
    return sum(counts)-max(counts)

def calcAlleleFreq(Allele,Reads):
    alleles =['A','C','G','T']
    Allele=alleles[Allele]
    tot=0
    if len(Reads.base)>0:
        for i in range(len(Reads.base)):
            if Reads.base[i]==Allele:
                tot+=1
    return tot

def calcGenoLogLike1(reads,site):
    alleles=['A','C','G','T']
    log_likes=[0.0,0.0,0.0,0.0,0.0]
    phredScale=33

    #cycle across all possible genotypes
    for j in range(len(alleles)):
        if j == 0:                
            for i in range(len(reads.base)):
            
         
                #get base probability from quality score
                bP = 10**((phredScale-ord(str(reads.base_quality[i])))/10)

                sublike=0.0
            
                if alleles[j]==reads.base[i]:
                    sublike += 1-(bP)
                else:
                    sublike += (bP/3)

                log_likes[j] += math.log(sublike)
                log_likes[4] += math.log(bP/3)
        else:
            for i in range(len(reads.base)):
                #get base probability from quality score
                bP = 10**((phredScale-ord(str(reads.base_quality[i])))/10)
                sublike=0.0
                if alleles[j]==reads.base[i]:
                    sublike += 1-(bP)
                else:
                    sublike += (bP/3)
                log_likes[j] += math.log(sublike)
    return log_likes

def exp_or_zero(x):
    if(x==0):
        x=0.0
    else:
        x=math.exp(x)
    return(x)

def log_or_zero(x):
    if(x==0):
        x=0.0
    else:
        x=math.log(x)
    return(x)

def delta_to_ploidy(delta_prob):
    output=0
    check=0
    for i in range(len(delta_prob)-1): 
        if math.exp(delta_prob[i])-math.exp(delta_prob[i+1])<0: # look for a case of where the jump in delta_prob has increased 
            output=i
            check=1
            break
    if check==0:
        output=1
    else:
        output=output
    return(output)

def dist(ploidies):
    ''' Function to return the distribution of ploidies predicted from the inputted array'''
    number = len(ploidies)
    dist = [ploidies.count(i)/number for i in range(1,9)]
    return dist


def filter(reads,min_quality_score):
    phredScale=33
    bases=""
    qualities=""
    for i in range(0,len(reads.base)):
        if ord(reads.base_quality[i])-phredScale>min_quality_score:
            bases+=reads.base[i]
            qualities+=reads.base_quality[i]
    return(bases,qualities)


# combinations with replacements, edited function from itertools
def combinations_with_rep(iterable, r):
    pool = list(iterable)
    n = len(pool)
    if not n and r:
        return
    indices = [0] * r
    yield list(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != n - 1:
                break
        else:
            return
        indices[i:] = [indices[i] + 1] * (r - i)
        yield list(pool[i] for i in indices)

# calculate genotype likelihoods (in ln format) in case of all given Nploids for Major and Minor
def calcGenoLogLikeN_MajorMinor(N,read,site,major,minor):
    alleles = ['A','C','G','T']
    log_likes=[0.0]*(N+1)
    it = -1
    phredScale=33
    mm = [major,minor]
    mmList = list(combinations_with_rep(mm,N)) # List of major minor combinations
    # cycle across all possible genotypes
    readLen = len(read.base)
    for subList in mmList:
        it += 1
        for i in range(readLen):
            bP = 10**((phredScale-ord(str(read.base_quality[i])))/10)
            sublike = 0.0
            for item in subList:
                if alleles[item] == read.base[i]:
                    sublike += (1-bP)/N
                else:
                    sublike += (bP/3)/N
            log_likes[it] += math.log(sublike)
    return(log_likes)
