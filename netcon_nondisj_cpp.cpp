// netcon_nondisj_cpp v2.01_arXiv
// by Robert N. C. Pfeifer, 2014
// Contact: rpfeifer.public@gmail.com, robert.pfeifer@mq.edu.au

// Note: This code uses the nomenclature of Appendix G of arXiv:1304.6112v4 [cond-mat.str-el] (2014) [or equivalently, Appendix F of Phys. Rev. E 90, 033315 (2014)], where a composite tensor which may be contracted with the result of an outer product is denoted X.
// Note 2: This code contains minor bugfixes when compared with v2.00 which was released with Phys. Rev. E 90, 033315 (2014) - see comments in netcon.m.

#include "mex.h"

class objectData
{
public:
    objectData();
    ~objectData();
    
    bool *legLinks;
    bool *tensorflags; // Length is always numtensors
    mwIndex sequencelen;
    mwIndex *sequence;
    mwIndex costlen;
    double *cost;
    bool isOP;
    double *maxdim; // Length 2 if costType==2, length 1 if costType==1
    double *allIn;
    bool isnew;
};

objectData::objectData()
{
    legLinks = NULL;
    tensorflags = NULL;
    sequence = NULL;
    cost = NULL;
    maxdim = NULL;
    allIn = NULL;
}
    
objectData::~objectData()
{
    if (legLinks!=NULL) {delete[] legLinks; legLinks = NULL;}
    if (tensorflags!=NULL) {delete[] tensorflags; tensorflags = NULL;}
    if (sequence!=NULL) {delete[] sequence; sequence = NULL;}
    if (cost!=NULL) {delete[] cost; cost = NULL;}
    if (maxdim!=NULL) {delete[] maxdim; maxdim = NULL;}
    if (allIn!=NULL) {delete[] allIn; allIn = NULL;}
}

// Objects are indexed in two different ways:
// (objectList*)objects[n] is a list of all objects made up of exactly n tensors.
// (ObjectTree*)root points to a tree structure used to quickly locate an object containing a specific set of tensors.


class objectList
{
public:
    objectList();
    ~objectList();
    
    objectList* next;
    objectData* object;
};

objectList::objectList()
{
    next = NULL;
    object = NULL;
}

objectList::~objectList()
{
    // Iterative deletion of downstream list (not recursive - lists can be hundreds of thousands of entries long, and recursion can overflow the stack)
    // Entire list is deleted when root node is deleted.
    while (next!=NULL) {
        objectList *t = next;
        next = next->next;
        t->next = NULL;
        delete t;
    }
    // Do not delete objects - this is done when deleting the object tree
    object = NULL;
}


class objectTree
{
public:
    objectTree(const mwIndex numtensors, const mwIndex offset);
    ~objectTree();
    inline objectTree* getbranch(const mwIndex branch);
    inline void setbranch(const mwIndex branch, objectTree *target);
    
    mwIndex numentries;
    objectData* object;
    mwIndex offset;
    
private:
    objectTree** branches;
};

objectTree::objectTree(const mwIndex numtensors, const mwIndex o)
{
    offset = o;
    numentries = numtensors-offset;
    branches = new objectTree*[numentries]; // OK
    for (mwIndex x=0;x<numentries;x++) branches[x] = NULL;
    object = NULL;
}

objectTree::~objectTree()
{
    // Recursive deletion of object tree (acceptable: max recursion depth = numtensors+1)
    // Entire tree is deleted when root node is deleted.
    for (mwIndex x=0;x<numentries;x++) if (branches[x]!=NULL) {delete branches[x]; branches[x] = NULL;}
    delete[] branches;
    branches = NULL;
    if (object!=NULL) {delete object; object=NULL;}
}

inline objectTree* objectTree::getbranch(const mwIndex branch) {return branches[branch-offset];}
inline void objectTree::setbranch(const mwIndex branch, objectTree *target) {branches[branch-offset] = target;}


class tensorXlist
{
public:
    tensorXlist();
    ~tensorXlist();
    
    tensorXlist *next;
    tensorXlist *prev;
    bool *legs;
    double *allIn; // For a tensor whose construction is consistent with Fig.5(c), allIn is the dimension of index c, corresponding to |E|.
    mwIndex flag;
};

tensorXlist::tensorXlist()
{
    next = NULL;
    prev = NULL;
    legs = NULL;
    allIn = NULL;
    flag = 2;
}

tensorXlist::~tensorXlist()
{
    // Only this node is deleted.
    next = NULL;
    prev = NULL;
    if (legs!=NULL) {
        delete[] legs;
        legs = NULL;
    }
    if (allIn!=NULL) {
        delete[] allIn;
        allIn = NULL;
    }
}


inline bool anyand(const bool *in1,const bool *in2,const mwIndex len) {for (mwIndex x=0;x<len;x++) if (in1[x] && in2[x]) return true; return false;}
inline bool any(const bool *in,const mwIndex len) {for (mwIndex x=0;x<len;x++) if (in[x]) return true; return false;}
inline bool notany(const bool *in,const mwIndex len) {for (mwIndex x=0;x<len;x++) if (in[x]) return false; return true;}
inline bool allAofB(const bool *in1,const bool *in2,const mwIndex len) {for (mwIndex x=0;x<len;x++) if (in2[x]) if (!in1[x]) return false; return true;}
inline bool* do_and(const bool *in1,const bool *in2,const mwIndex len) {bool *rtn = new bool[len]; for (mwIndex x=0;x<len;x++) rtn[x] = in1[x] && in2[x]; return rtn;} // *
inline bool* do_andnot(const bool *in1,const bool *in2,const mwIndex len) {bool *rtn = new bool[len]; for (mwIndex x=0;x<len;x++) rtn[x] = in1[x] && (~in2[x]); return rtn;} // *
inline bool* do_xor(const bool *in1,const bool *in2,const mwIndex len) {bool *rtn = new bool[len]; for (mwIndex x=0;x<len;x++) rtn[x] = (in1[x] && !in2[x]) || (!in1[x] && in2[x]); return rtn;} // *
inline bool* do_or(const bool *in1,const bool *in2,const mwIndex len) {bool *rtn = new bool[len]; for (mwIndex x=0;x<len;x++) rtn[x] = (in1[x] || in2[x]); return rtn;} // *
inline bool equalflags(const bool *in1,const bool *in2,const mwIndex len) {for (mwIndex x=0;x<len;x++) if (in1[x]!=in2[x]) return false; return true;}
inline void pause(const double pauselen) {mxArray *pause_mx = mxCreateDoubleScalar(pauselen); mxArray *prhs[] = {pause_mx}; mexCallMATLAB(0,NULL,1,prhs,"pause"); mxDestroyArray(pause_mx); pause_mx = NULL;}

inline double* getprodlegdims(const bool *freelegs, const double *legCosts, const mwSize numleglabels, const mwIndex costType); // *
inline double* dimSquared(const double *dim, const mwIndex costType); // *
inline bool isgreaterthan_sd(const double *dim1, const double *dim2, const mwIndex costType);
inline double* getbuildcost(const bool *freelegs, const bool *commonlegs, const double *legCosts, const mwIndex costType, const double oldmuCap, const double muCap, bool isnew,
        const double *cost1, const double *cost2, const mwIndex cost1len, const mwIndex cost2len, const mwIndex numleglabels, double &rtnmuCap, bool &isOK, mwIndex &newCostLen); // *
inline bool islessthan(const double *cost1,const double *cost2,const mwIndex len1,const mwIndex len2,const mwIndex costType);
inline void displaycostandsequence(const objectData *ptr, const double costtype, const mxArray* tracedindices, const mxArray *posindices);
inline void addToTensorXlist(tensorXlist *&Xlistroot, const bool *freelegs, const mwIndex numleglabels, const double *newallin, const mwIndex costType, const bool isnew, const bool oldMayHaveEntry);
inline void updateTensorXlist(tensorXlist *&Xlistroot, const bool *freelegs, const mwIndex numleglabels, const double *newallin, const mwIndex costType);
inline void removeFromTensorXlist(tensorXlist *&Xlistroot, const bool *freelegs, const mwIndex numleglabels);
inline void mergeTensorXlist(tensorXlist *&Xlistroot, const mwIndex costType, const mwIndex numleglabels);
inline bool allundertwo(tensorXlist *list);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // Main routine
    // [sequence cost] = netcon_nondisj_cpp(legLinks,legCosts,verbosity,costType,muCap,allowOPs,posindices,tracedindices)
    
    // Extract supplied data
    // =====================
    const mxArray *legLinks = prhs[0];
    const double *legCosts = mxGetPr(prhs[1]);
    const double verbosity = mxGetScalar(prhs[2]);
    const double costType = mxGetScalar(prhs[3]);
    double muCap = mxGetScalar(prhs[4]);
    const bool allowOPs = (mxGetScalar(prhs[5])==1);
    const mxArray *posindices = prhs[6];
    const mxArray *tracedindices = prhs[7];
    
    // Prepare useful vars
    // ===================
    const mwSize *legCostsdims = mxGetDimensions(prhs[1]);
    const mwSize numleglabels = legCostsdims[0];
    const mwIndex numtensors = mxGetNumberOfElements(legLinks);
    double minlegcost;
    if (muCap > 1e300) mexPrintf("Warning from netcon_nondisj_cpp: muCap very large - risk of numerical overflow.\nRecommended action: Reduce muCap and try again.\n\n");
    if (costType==1) {
        minlegcost = legCosts[0];
        for (mwIndex a=1;a<numleglabels;a++) if (legCosts[a]<minlegcost) minlegcost = legCosts[a];
        if (minlegcost<2) minlegcost = 2;
        if (muCap<minlegcost) muCap=minlegcost;
    }
    double oldmuCap = 0;
    double newmuCap = 1e308;
    bool done = false;
    
    // Create objects index
    // ====================
    objectTree *root = new objectTree(numtensors,0); // OK
    objectList **objects = new objectList*[numtensors]; // OK
    objectList *newObjects = NULL;
    for (mwIndex a=0;a<numtensors;a++) objects[a] = NULL;
    
    // ### Create linked lists used in enforcing Sec. II.B.2.c (structure of tensor X contractable with an outer product)
    tensorXlist *Xlistroot = new tensorXlist; // OK
    tensorXlist *Xlistptr = Xlistroot;
    // ### End of creating linked lists

    for (mwIndex a=0;a<numtensors;a++) {
        root->setbranch(a,new objectTree(numtensors,0)); // OK
        root->getbranch(a)->object = new objectData; // OK
        int *legLinksdata = (int*)mxGetPr(mxGetCell(legLinks,a));
        root->getbranch(a)->object->legLinks = new bool[numleglabels]; // OK
        for (mwIndex b=0;b<numleglabels;b++) root->getbranch(a)->object->legLinks[b] = false;
        for (mwIndex b=0;b<mxGetNumberOfElements(mxGetCell(legLinks,a));b++)
            root->getbranch(a)->object->legLinks[legLinksdata[b]-1] = true; // Leg labels are 1-based
        legLinksdata = NULL;

        // ### Set up initial linked list data used in enforcing Sec. II.B.2.c
        if (a!=0) {
            Xlistptr->next = new tensorXlist; // OK
            Xlistptr->next->prev = Xlistptr;
            Xlistptr = Xlistptr->next;
        }
        Xlistptr->legs = new bool[numleglabels]; // OK
        for (mwIndex b=0;b<numleglabels;b++) Xlistptr->legs[b] = root->getbranch(a)->object->legLinks[b];
        Xlistptr->allIn = new double[mwIndex(costType)]; // OK
        for (mwIndex b=0;b<costType;b++) Xlistptr->allIn[b] = 1e308;
        Xlistptr->flag = 0;
        // ### End of setting up initial linked list data

        root->getbranch(a)->object->tensorflags = new bool[numtensors]; // OK
        for (mwIndex b=0;b<numtensors;b++) root->getbranch(a)->object->tensorflags[b] = (b==a);
        root->getbranch(a)->object->sequencelen = 0;
        root->getbranch(a)->object->costlen = 1;
        root->getbranch(a)->object->cost = new double[1]; // OK
        root->getbranch(a)->object->cost[0] = 0;
        root->getbranch(a)->object->isOP = false;
        root->getbranch(a)->object->allIn = new double[mwIndex(costType)]; // OK
        root->getbranch(a)->object->isnew = true;
        for (mwIndex b=0;b<costType;b++) root->getbranch(a)->object->allIn[b] = 0;

        objectList *t = new objectList; // OK
        t->next = objects[0];
        t->object = root->getbranch(a)->object;
        objects[0] = t;
        t = new objectList; // OK
        t->next = newObjects;
        t->object = root->getbranch(a)->object;
        newObjects = t;
    }
    
    // Main loop
    // =========
    while (!done) {
        if (verbosity>0 && oldmuCap!=muCap) {
            if (costType==1) {
                mexPrintf("Looking for solutions with maximum cost of %g\n",muCap);
                pause(0.01);
            } else {
                mexPrintf("Looking for solutions of cost O(X^%g)\n",muCap);
                pause(0.01);
            }
        }
        for (mwIndex numInObjects=2;numInObjects<=numtensors;numInObjects++) { // All numIn values are 1-based, so note -1 when used as an index
            if (verbosity > 2) mexPrintf("Pairwise contractions (AB) involving %d fundamental tensors:\n",numInObjects);
            for (mwIndex numInPieceOne=1;numInPieceOne<=numInObjects/2;numInPieceOne++) { // mwIndex is integer so effectively rounds down
                mwIndex numInPieceTwo = numInObjects-numInPieceOne;
                if (verbosity > 2) {
                    mexPrintf("A contains %d, B contains %d.\n",numInPieceOne,numInPieceTwo);
                    pause(0.01);
                }
                if (objects[numInPieceOne-1]!=NULL && objects[numInPieceTwo-1]!=NULL) {
                    // Iterate over pairings: Iterate over object 1
                    for (objectList *a = objects[numInPieceOne-1];a!=NULL;a=a->next) {
                        // Get data pointer for object 1
                        objectData *objectPtr1 = a->object;
                        // obj1data: legflags, tensorflags, sequence, cost(ToBuild), isOP, maxdim(encounteredDoingOPs), allIn
                        const bool isnew1 = objectPtr1->isnew;
                        
                        // Iterate over pairings: Iterate over object 2
                        for (objectList *b = (numInPieceOne==numInPieceTwo) ? a->next : objects[numInPieceTwo-1];b!=NULL;b=b->next) {
                            // Check object 1 and object 2 don't share any common tensors (which would then appear twice in the resulting network)
                            // Get data pointer for object 2
                            objectData *objectPtr2 = b->object;
                            if (!anyand(objectPtr1->tensorflags,objectPtr2->tensorflags,numtensors)) {
                                // obj2data: legflags, tensorflags, sequence, cost(ToBuild), isOP, maxdim(encounteredDoingOPs), allIn
                                const bool isnew2 = objectPtr2->isnew;
                                
                                const bool *legs1 = objectPtr1->legLinks;
                                const bool *legs2 = objectPtr2->legLinks;
                                bool *freelegs = do_xor(legs1,legs2,numleglabels); // OK
                                bool *freelegs1 = do_and(legs1,freelegs,numleglabels); // OK
                                bool *freelegs2 = do_and(legs2,freelegs,numleglabels); // OK
                                bool *commonlegs = do_and(legs1,legs2,numleglabels); // OK
                                bool commonlegsflag = any(commonlegs,numleglabels);
                                
                                // Exclude outer products if allowOPs not set
                                bool isOK = (allowOPs || commonlegsflag);
                                int thisTensorXflag = -1; // If contracting an outer product and a tensor X, this records the flag associated with that tensor X. -1 is a recogniseable dummy value.
                                
                                // ### Enforce Sec. II.B.2.b,c,d (only perform outer product if there is a known tensor X with appropriate structure; only contract resulting object in another outer product or with an appropriate tensor X; enforce index dimension constraints)
                                if (isOK && !commonlegsflag) {
                                    // It's an outer product. Check if a suitable tensor X exists yet to contract with this outer product [Fig.5(c) & Eq.(25)].
                                    mwIndex flagStop;
                                    if (oldmuCap==muCap) {
                                        // Pass for new X's
                                        if (isnew1 || isnew2) {
                                            // Using a new object - allowed to contract with old and new X's
                                            Xlistptr = Xlistroot;
                                            flagStop = 1;
                                        }
                                        else {
                                            // Made from old objects - only allowed to contract with new X's
                                            Xlistptr = Xlistroot;
                                            bool done = false;
                                            while (!done) {
                                                if (Xlistptr==NULL) done = true;
                                                else {
                                                    if (Xlistptr->flag==0) Xlistptr = Xlistptr->next;
                                                    else done = true;
                                                }
                                            }
                                            flagStop = 1;
                                        }
                                    }
                                    else {
                                        // Old X's only on this pass
                                        Xlistptr = Xlistroot;
                                        flagStop = 0;
                                    }
                                    double *tdim = NULL;
                                    double *tdim1 = NULL;
                                    double *tdim2 = NULL;
                                    while (Xlistptr!=NULL) {
                                        if (Xlistptr->flag > flagStop) Xlistptr = NULL;
                                        else {
                                            if (allAofB(Xlistptr->legs,freelegs,numleglabels)) {
                                                // IIB2c: xi_C > xi_A (25)
                                                if (tdim==NULL) tdim = getprodlegdims(freelegs,legCosts,numleglabels,costType); // OK
                                                if (isgreaterthan_sd(Xlistptr->allIn,tdim,costType)) {
                                                    // IIB2b: xi_C > xi_D && xi_C > xi_E: (16)
                                                    if (tdim1==NULL) tdim1 = getprodlegdims(freelegs1,legCosts,numleglabels,costType); // OK
                                                    bool *Xfreelegs = do_andnot(Xlistptr->legs,freelegs,numleglabels); // OK
                                                    double *tdimX = getprodlegdims(Xfreelegs,legCosts,numleglabels,costType); // OK
                                                    if (isgreaterthan_sd(tdimX,tdim1,costType)) {
                                                        if (tdim2==NULL) tdim2 = getprodlegdims(freelegs2,legCosts,numleglabels,costType); // OK
                                                        if (isgreaterthan_sd(tdimX,tdim2,costType)) {
                                                            thisTensorXflag = Xlistptr->flag;
                                                            delete[] Xfreelegs; Xfreelegs = NULL;
                                                            delete[] tdimX; tdimX = NULL;
                                                            break;
                                                        }
                                                    }
                                                    delete[] Xfreelegs; Xfreelegs = NULL;
                                                    delete[] tdimX; tdimX = NULL;
                                                }
                                            }
                                            Xlistptr = Xlistptr->next;
                                        }
                                    }
                                    if (tdim2!=NULL) {delete[] tdim2; tdim2 = NULL;}
                                    if (tdim1!=NULL) {delete[] tdim1; tdim1 = NULL;}
                                    if (tdim!=NULL) {delete[] tdim; tdim = NULL;}
                                    isOK = (thisTensorXflag!=-1); // Outer products only OK if a corresponding tensor X has been found
                                }
                                
                                // If either constituent is the result of an outer product, check that it is being contracted with an appropriate tensor
                                // [either this is a contraction over all indices, or this is an outer product with a tensor of larger total dimension than
                                // either constituent of the previous outer product, and Eqs. (16), (25), and (27) are satisfied].
                                if (isOK && (objectPtr1->isOP || objectPtr2->isOP)) {
                                    // Post-OP. This contraction only allowed if it is also an outer product, or if only one object is an outer product, it involves all indices on that tensor, and the other object satisfies the relevant conditions.
                                    isOK = (!objectPtr1->isOP) || (!objectPtr2->isOP) || (!commonlegsflag); // If contracting over common indices, only one object may be an outer product
                                    if (isOK) {
                                        if (commonlegsflag) {
                                            // This contraction is not itself an outer product
                                            double *freelegsdim = NULL;
                                            // Conditions on outer product object:
                                            if (objectPtr1->isOP) { // Object 1 is an outer product
                                                // Check all-indices condition:
                                                if (equalflags(commonlegs,legs1,numleglabels)) {
                                                    // Check free legs on contracting tensor are larger than summing legs going to each component of outer product [Eq. (16)]
                                                    freelegsdim = getprodlegdims(freelegs,legCosts,numleglabels,costType); // OK
                                                    isOK = isgreaterthan_sd(freelegsdim,objectPtr1->maxdim,costType); // IIB2b: xi_C > xi_D, xi_C > xi_E (16)
                                                }
                                                else isOK = false;
                                            }
                                            else { // Object 2 is an outer product
                                                // Check all-indices condition:
                                                if (equalflags(commonlegs,legs2,numleglabels)) {
                                                    // Check free legs on contracting tensor are larger than summing legs going to each component of outer product [Eq. (16)]
                                                    freelegsdim = getprodlegdims(freelegs,legCosts,numleglabels,costType); // OK
                                                    isOK = isgreaterthan_sd(freelegsdim,objectPtr2->maxdim,costType); // IIB2b: xi_C > xi_D, xi_C > xi_E (16)
                                                }
                                                else isOK = false;
                                            }
                                            // Conditions on X: Ensure X is fundamental or acceptably-constructed (note: structure is checked by requiring
                                            // non-zero value of allIn)
                                            if (isOK) {
                                                if (objectPtr1->isOP) {
                                                    // Tensor 2 is X
                                                    if (numInPieceTwo > 1) {
                                                        // Tensor 2 is not fundamental
                                                        // Check tensor 2 is constructed in an acceptable fashion [Fig. 5(c) and Eqs. (25) and (26)]
                                                        isOK = isgreaterthan_sd(objectPtr2->allIn,freelegsdim,costType); // IIB2c: xi_C > xi_D (26)
                                                        if (isOK) {
                                                            delete[] freelegsdim;
                                                            freelegsdim = getprodlegdims(freelegs1,legCosts,numleglabels,costType); // OK
                                                            isOK = isgreaterthan_sd(objectPtr2->allIn,freelegsdim,costType); // IIB2c: xi_C > xi_A (25)
                                                        }
                                                    }
                                                }
                                                else {
                                                    // Tensor 1 is X
                                                    if (numInPieceOne > 1) {
                                                        // Tensor 1 is not fundamental
                                                        // Check tensor 1 is constructed in an acceptable fashion [Fig. 5(c) and Eqs. (25) and (26)]
                                                        isOK = isgreaterthan_sd(objectPtr1->allIn,freelegsdim,costType); // IIB2c: xi_C > xi_D (26)
                                                        if (isOK) {
                                                            delete[] freelegsdim;
                                                            freelegsdim = getprodlegdims(freelegs2,legCosts,numleglabels,costType); // OK
                                                            isOK = isgreaterthan_sd(objectPtr1->allIn,freelegsdim,costType); // IIB2c: xi_C > xi_A (25)
                                                        }
                                                    }
                                                }
                                            }
                                            if (freelegsdim!=NULL) delete[] freelegsdim;
                                        }
                                        else {
                                            // This contraction is an outer product. If either constituent is an outer product, check that both tensors
                                            // within that object are not larger than the third tensor with which they are now being contracted.
                                            if (objectPtr1->isOP) {
                                                double *freelegsdim = getprodlegdims(freelegs2,legCosts,numleglabels,costType); // OK
                                                isOK = !isgreaterthan_sd(objectPtr1->maxdim,freelegsdim,costType); // IIB2b: xi_C >= xi_A, xi_C >= xi_B (20)
                                                delete[] freelegsdim;
                                            }
                                            if (isOK) {
                                                if (objectPtr2->isOP) {
                                                    double *freelegsdim = getprodlegdims(freelegs1,legCosts,numleglabels,costType); // OK
                                                    isOK = !isgreaterthan_sd(objectPtr2->maxdim,freelegsdim,costType); // IIB2b: xi_C >= xi_A, xi_C >= xi_B (20)
                                                    delete[] freelegsdim;
                                                }
                                            }
                                        }
                                    }
                                }
                                // ### End of enforcing Sec. II.B.2.b,c,d  (only perform outer product if there is a known tensor X with appropriate structure; only contract resulting object in another outer product or with an appropriate tensor X; enforce index dimension constraints)
                                
                                // If contraction is not prohibited, check cost is acceptable (<=muCap and, if not involving new objects, >oldmuCap)
                                double *newCost = NULL;
                                mwIndex newCostLen;
                                if (isOK) {
                                    // If constructing an outer product which may contract with a new X, do not exclude on basis of low cost: Hence
                                    // isnew1||isnew2||thisTensorXflag>0
                                    double rtnmuCap;
                                    newCost = getbuildcost(freelegs,commonlegs,legCosts,costType,oldmuCap,muCap,isnew1||isnew2||(thisTensorXflag>0),objectPtr1->cost,objectPtr2->cost,objectPtr1->costlen,objectPtr2->costlen,numleglabels,rtnmuCap,isOK,newCostLen); // OK
                                    if (!isOK) if (rtnmuCap < newmuCap && rtnmuCap!=0) newmuCap = rtnmuCap;
                                }

                                // If cost is OK, compare with previous best known cost for constructing this object
                                bool *tensorsInNew = NULL;
                                objectData *objptr = NULL;
                                bool isnew;
                                double *newallin = NULL;
                                if (isOK) {
                                    // Get involved tensors
                                    tensorsInNew = do_or(objectPtr1->tensorflags,objectPtr2->tensorflags,numtensors); // OK
                                    // Find if previously constructed
                                    objectTree *treeptr = root;
                                    bool isnew = false;
                                    for (mwIndex x=0;x<numtensors;x++) if (tensorsInNew[x]) {
                                        if (treeptr->getbranch(x)==NULL) {
                                            isnew = true;
                                            treeptr->setbranch(x,new objectTree(numtensors,x+1)); // OK
                                        }
                                        treeptr = treeptr->getbranch(x);
                                    }
                                    if (treeptr->object==NULL) isnew = true;
                                    if (isnew) {
                                        // Create space for object data
                                        treeptr->object = new objectData; // OK
                                        // Add to list of objects of this size
                                        objectList *t = new objectList; // OK
                                        t->object = treeptr->object;
                                        t->next = objects[numInObjects-1];
                                        objects[numInObjects-1] = t;
                                    }
                                    objptr = treeptr->object;
                                    if (!isnew) {
                                        // Compare new cost with best-so-far cost for construction of this object
                                        isOK = islessthan(newCost,objptr->cost,newCostLen,objptr->costlen,costType);
                                    }

                                    // ### If appropriate, update tensorXlist (list of tensors which can be contracted with objects created by outer product)
                                    if (allowOPs) {
                                        if (isOK) {
                                            // New tensor or new best cost
                                            bool E_is_1 = (notany(freelegs1,numleglabels) && any(freelegs2,numleglabels));
                                            bool E_is_2 = (notany(freelegs2,numleglabels) && any(freelegs1,numleglabels));
                                            if (E_is_1 || E_is_2) {
                                                // New best sequence consistent with Fig.5(c).
                                                // Determine the value of allIn, which corresponds to xi_C. (This is used in determining valid tensors X to contract with outer products).
                                                if (E_is_1) newallin = getprodlegdims(legs1,legCosts,numleglabels,costType); // OK
                                                else newallin = getprodlegdims(legs2,legCosts,numleglabels,costType); // OK
                                                // Add to tensor X list for outer products (or if already there, update the value of allIn):
                                                double *dimSq = dimSquared(newallin,costType); // OK
                                                double *freelegsdim = getprodlegdims(freelegs,legCosts,numleglabels,costType); // OK
                                                if (isnew) {
                                                    if (isgreaterthan_sd(dimSq,freelegsdim,costType)) { // Enforce Eq.(27)
                                                        addToTensorXlist(Xlistroot,freelegs,numleglabels,newallin,costType,true,false); // Permitted to act as a tensor X: Add to list
                                                    }
                                                    else {
                                                        // Set allIn to zero. While it is not actually zero, and a lower value may be found later, this is
                                                        // irrelevant as neither this value nor any lower one are compatible with the constraints of
                                                        // Eq.(27). This tensor never acts as a tensor X.
                                                        for (mwIndex x=0;x<costType;x++) newallin[x] = 0;
                                                    }
                                                }
                                                else {
                                                    bool oldMayHaveEntry;
                                                    if (costType==1) oldMayHaveEntry = (objptr->allIn[0]!=0);
                                                    else oldMayHaveEntry = (objptr->allIn[0]!=0) || (objptr->allIn[1]!=0);
                                                    if (isgreaterthan_sd(dimSq,freelegsdim,costType)) { // Enforce Eq.(27)
                                                        addToTensorXlist(Xlistroot,freelegs,numleglabels,newallin,costType,false,oldMayHaveEntry); // Permitted to act as a tensor X: Add to list
                                                    }
                                                    else {
                                                        // Set allIn to zero. While it is not actually zero, and a lower value may be found later, this is
                                                        // irrelevant as neither this value nor any lower one are compatible with the constraints of
                                                        // Eq.(27). This tensor never acts as a tensor X.
                                                        for (mwIndex x=0;x<costType;x++) newallin[x] = 0;
                                                        // Also remove from tensor X list if present.
                                                        if (oldMayHaveEntry) removeFromTensorXlist(Xlistroot,freelegs,numleglabels);
                                                    }
                                                }
                                                delete[] freelegsdim;
                                                delete[] dimSq;
                                            }
                                            else {
                                                // This tensor is not an eligible tensor X for an outer product: Store a dummy value in allIn to indicate this
                                                newallin = new double[mwIndex(costType)]; // OK
                                                for (mwIndex x=0;x<costType;x++) newallin[x] = 0;
                                                // Best cost and not consistent with Fig.5(c): Ensure does not appear in tensorXlist. Active removal only
                                                // required if object is not new, and previous best sequence was consistent with Fig.5(c), so allIn is not a
                                                // dummy on the old entry.
                                                if (!isnew) {
                                                    bool oldMayHaveEntry;
                                                    if (costType==1) oldMayHaveEntry = (objptr->allIn[0]!=0);
                                                    else oldMayHaveEntry = (objptr->allIn[0]!=0) || (objptr->allIn[1]!=0);
                                                    if (oldMayHaveEntry) removeFromTensorXlist(Xlistroot,freelegs,numleglabels);
                                                }
                                            }
                                        }
                                        else {
                                            bool equalCosts = objptr->costlen==newCostLen;
                                            if (equalCosts) {
                                                for (mwIndex x=0;x<newCostLen;x++) if (objptr->cost[x]!=newCost[x]) {
                                                    equalCosts = false;
                                                    break;
                                                }
                                            }
                                            if (equalCosts) {
                                                // Equal-best cost to a known sequence for the same tensor
                                                bool oldMayHaveEntry;
                                                if (costType==1) oldMayHaveEntry = (objptr->allIn[0]!=0);
                                                else oldMayHaveEntry = (objptr->allIn[0]!=0) || (objptr->allIn[1]!=0);
                                                if (oldMayHaveEntry) {
                                                    // Previous best sequence was consistent with Fig.5(c) so tensor may appear in the provisional environments list
                                                    bool E_is_1 = (notany(freelegs1,numleglabels) && any(freelegs2,numleglabels));
                                                    bool E_is_2 = (notany(freelegs2,numleglabels) && any(freelegs1,numleglabels));
                                                    if (E_is_1 || E_is_2) {
                                                        // Determine the value of allIn, which corresponds to xi_C in Fig.5(c).
                                                        if (E_is_1) newallin = getprodlegdims(legs1,legCosts,numleglabels,costType); // OK
                                                        else newallin = getprodlegdims(legs2,legCosts,numleglabels,costType); // OK
                                                        // If smaller than previous value, update the value of allIn:
                                                        if (isgreaterthan_sd(objptr->allIn,newallin,costType)) {
                                                            double *dimSq = dimSquared(newallin,costType); // OK
                                                            double *freelegsdim = getprodlegdims(freelegs,legCosts,numleglabels,costType); // OK
                                                            if (isgreaterthan_sd(dimSq,freelegsdim,costType)) { // Enforce Eq.(27)
                                                                updateTensorXlist(Xlistroot,freelegs,numleglabels,newallin,costType); // Update entry in tensor X list
                                                                double *t = objptr->allIn; // Update minimum value of allIn in object data
                                                                objptr->allIn = newallin;
                                                                newallin = t;
                                                            }
                                                            else {
                                                                removeFromTensorXlist(Xlistroot,freelegs,numleglabels); // Remove entry from tensor X list
                                                                for (mwIndex x=0;x<costType;x++) objptr->allIn[x] = 0; // Zero minimum value of allIn in object data [never acts as a tensor X, by Eq.(27)]
                                                            }
                                                            delete[] dimSq;
                                                            delete[] freelegsdim;
                                                        }
                                                        delete[] newallin;
                                                        newallin = NULL;
                                                    }
                                                    else {
                                                        // Found best-equal sequence not consistent with Fig.5(c)
                                                        removeFromTensorXlist(Xlistroot,freelegs,numleglabels);
                                                        for (mwIndex x=0;x<costType;x++) objptr->allIn[x] = 0; // Zero minimum value of allIn in object data (never acts as a tensor X)
                                                    }
                                                }
                                                //else: There already exists a best-known-cost sequence for the tensor which is not consistent with Fig.5(c), and isOK=false. Tensor does not appear in tensorXlist. No need to assign allIn. Now returning to start of main loop.
                                            }
                                            //else: Sequence is not capable of updating tensorXlist (not better cost, not equal cost). Also, isOK=false. No need to assign allIn. Now returning to start of main loop.
                                        }
                                    }
                                    //else: Not doing outer products. Leave NULL as a dummy value in allIn, which is never used.
                                    // ### Done updating tensorXlist (list of tensors which can be contracted with objects created by outer product)
                                }

                                
                                if (isOK) {
                                    // Either no previous construction, or this one is better
                                    // Update object data with this construction
                                    
                                    // Ensure that when compositing sequences, the outer product object (if there is one) goes second:
                                    objectData *seq1obj, *seq2obj;
                                    if (objectPtr1->isOP) {
                                        seq1obj = objectPtr2;
                                        seq2obj = objectPtr1;
                                    }
                                    else {
                                        seq1obj = objectPtr1;
                                        seq2obj = objectPtr2;
                                    }
                                    
                                    // Record object
                                    if (!commonlegsflag) {
                                        // ### This construction is an outer product. Note dimension of larger of the two participating tensors in newmaxdim. (This is used in enforcing index-dimension-related constraints.)
                                        double *newmaxdim = NULL;
                                        newmaxdim = getprodlegdims(legs1,legCosts,numleglabels,costType); // OK
                                        double *newmaxdim2 = getprodlegdims(legs2,legCosts,numleglabels,costType); // OK
                                        if (isgreaterthan_sd(newmaxdim2,newmaxdim,costType)) {
                                            delete[] newmaxdim;
                                            newmaxdim = newmaxdim2;
                                            newmaxdim2 = NULL;
                                        }
                                        else {
                                            delete[] newmaxdim2;
                                            newmaxdim2 = NULL;
                                        }
                                        if (objptr->maxdim!=NULL) delete[] objptr->maxdim;
                                        objptr->maxdim = newmaxdim;
                                        newmaxdim = NULL;
                                        // ### End recording dimension of larger of the two participating tensors in newmaxdim.
                                        if (objptr->sequence!=NULL) delete[] objptr->sequence;
                                        objptr->sequencelen = objectPtr1->sequencelen + objectPtr2->sequencelen + 1;
                                        objptr->sequence = new mwIndex[objptr->sequencelen]; // OK
                                        mwIndex ptr=0;
                                        for (mwIndex x=0;x<seq1obj->sequencelen;x++) objptr->sequence[ptr++] = seq1obj->sequence[x];
                                        for (mwIndex x=0;x<seq2obj->sequencelen;x++) objptr->sequence[ptr++] = seq2obj->sequence[x];
                                        objptr->sequence[ptr] = 0;
                                        objptr->isOP = true;
                                    }
                                    else {
                                        // This construction is not an outer product.
                                        mwIndex t = 0;
                                        for (mwIndex x=0;x<numleglabels;x++) if (commonlegs[x]) t++;
                                        if (objptr->sequence!=NULL) delete[] objptr->sequence;
                                        objptr->sequencelen = objectPtr1->sequencelen + objectPtr2->sequencelen + t;
                                        objptr->sequence = new mwIndex[objptr->sequencelen]; // OK
                                        t = 0;
                                        for (mwIndex x=0;x<seq1obj->sequencelen;x++) {objptr->sequence[t++] = seq1obj->sequence[x];}
                                        for (mwIndex x=0;x<seq2obj->sequencelen;x++) {objptr->sequence[t++] = seq2obj->sequence[x];}
                                        for (mwIndex x=0;x<numleglabels;x++) if (commonlegs[x]) {objptr->sequence[t++] = x+1;}
                                        objptr->isOP = false;
                                        // ### This construction is not an outer product. Therefore store a dummy value in maxdim. (For outer products, maxdim records the dimension of the larger participating tensor, to assist in enforcing index-dimension-related constraints.)
                                        if (objptr->maxdim!=NULL) delete[] objptr->maxdim;
                                        objptr->maxdim = NULL;
                                        // ### End storing dummy value in maxdim
                                    }
                                    if (objptr->legLinks!=NULL) delete[] objptr->legLinks;
                                    objptr->legLinks = freelegs;
                                    freelegs = NULL;
                                    if (objptr->tensorflags!=NULL) delete[] objptr->tensorflags;
                                    objptr->tensorflags = new bool[numtensors]; // OK
                                    for (mwIndex x=0;x<numtensors;x++) objptr->tensorflags[x] = tensorsInNew[x];
                                    if (objptr->cost!=NULL) delete[] objptr->cost;
                                    objptr->cost = newCost;
                                    objptr->costlen = newCostLen;
                                    newCost = NULL;
                                    // ### If this tensor has the structure of Fig.5(c) and so is capable of being contracted with an outer product object, |E| is recorded in newallin (otherwise this is a dummy value). Store this value.
                                    objptr->allIn = newallin;
                                    newallin = NULL;
                                    // ### Done storing value of |E| (if applicable).

                                    // Flag as new construction
                                    objptr->isnew = true;
                                    objectList *t = new objectList; // OK
                                    t->next = newObjects;
                                    t->object = objptr;
                                    newObjects = t;
                                    
                                    // If top level, display result
                                    if (numInObjects == numtensors) {
                                        // ### If a valid contraction sequence has been found, there is no need to perform any contraction sequence more expensive than this. Set muCap accordingly.
                                        if (costType==1) muCap = objptr->cost[0];
                                        else muCap = objptr->costlen-1;
                                        // ### Done setting muCap accordingly.
                                        if (verbosity > 1) displaycostandsequence(objptr,costType,tracedindices,posindices);
                                    }
                                }
                                if (freelegs!=NULL) {delete[] freelegs; freelegs = NULL;}
                                delete[] freelegs1; freelegs1 = NULL;
                                delete[] freelegs2; freelegs2 = NULL;
                                delete[] commonlegs; commonlegs = NULL;
                                if (newCost!=NULL) {delete[] newCost; newCost = NULL;}
                                if (tensorsInNew!=NULL) {delete[] tensorsInNew; tensorsInNew = NULL;}
                            }
                        }
                    }
                }
            }
        }
        
        // ### Check there are no new entries in the list of tensors which can be contracted with outer products.
        const bool allundertwoflag = allundertwo(Xlistroot);
        // ### Finished searching if an object has been constructed which contains all tensors, and no new outer products have been enabled on the last pass (allundertwoflag==true).
        done = objects[numtensors-1]!=NULL && (!allowOPs || allundertwoflag);
        
        if (!done) {
            if (allundertwoflag) {
                // ### All X tensors have been present for an entire pass, so all permitted outer products at this cost have already been constructed.
                // Increment muCap, update oldmuCap
                if (costType==1) if (newmuCap < muCap * minlegcost) newmuCap = muCap * minlegcost;
                oldmuCap = muCap;
                muCap = newmuCap;
                if (muCap > 1e300) mexPrintf("Warning from netcon_nondisj_cpp: muCap very large - risk of numerical overflow.\n\n");
                newmuCap = 1e308;
            }
            else {
                // ### New X tensors generated this pass (some tensor X flags==2). Do another pass with same cost limit, to construct newly-allowed objects (i.e. sequences of affordable cost including at least one new outer product).
                // ### This is achieved by updating oldmuCap only. Now only outer products and contractions involving newly-created tensors will satisfy mu_0 < mu <= muCap.
                oldmuCap = muCap;
            }
            // Clear all new object flags
            for (objectList *t = newObjects; t!=NULL; t=t->next) t->object->isnew = false;
            if (newObjects!=NULL) {delete newObjects; newObjects = NULL;}
            // ### Update tensor X flags (2 -> 1 -> 0):
            // ### 2: Newly created this pass becomes
            // ### 1: Created last pass; allow construction of cheap objects which contract with this, as they may have previously been excluded due to lack of a valid tensor X. 1 becomes...
            // ### 0: Old tensor X. Standard costing rules apply.
            // ### Delete redundant entries in tensorXlist (e.g. if A has a subset of the legs on B, and an equal or lower value of allIn (i.e. |E| in Fig.5(c)))
            mergeTensorXlist(Xlistroot,costType,numleglabels);
            // ### Done updating tensor X flags
        }
    }
    
    // Extract final result
    // ====================
    objectData *objectPtr = objects[numtensors-1]->object;
    mwSize *seqdims = new mwSize[2]; seqdims[0] = 1; seqdims[1] = objectPtr->sequencelen; // OK
    plhs[0] = mxCreateNumericArray(2,seqdims,mxINT32_CLASS,mxREAL);
    delete[] seqdims;
    seqdims = NULL;
    unsigned int *u32t_dbl = (unsigned int*)mxGetPr(plhs[0]);
    for (mwIndex x=0;x<objectPtr->sequencelen;x++) u32t_dbl[x] = objectPtr->sequence[x];
    plhs[1] = mxCreateDoubleMatrix(1,objectPtr->costlen,mxREAL);
    double *t_dbl;
    t_dbl = mxGetPr(plhs[1]);
    for (mwIndex x=0;x<objectPtr->costlen;x++) t_dbl[x] = objectPtr->cost[x];
    objectPtr = NULL;

    // Tidy up
    // =======
    while (Xlistroot!=NULL) {
        Xlistptr = Xlistroot;
        Xlistroot = Xlistroot->next;
        delete Xlistptr;
    }
    if (newObjects!=NULL) {delete newObjects; newObjects=NULL;}
    for (mwIndex x=0;x<numtensors;x++) {delete objects[x]; objects[x] = NULL;}
    delete[] objects;
    delete root;
}    
    



inline double* getprodlegdims(const bool *freelegs, const double *legCosts, const mwSize numleglabels, const mwIndex costType)
{
    double *dim;
    if (costType==1) {
        dim = new double[1]; // *
        dim[0] = 1;
        for (mwIndex a=0;a<numleglabels;a++) if (freelegs[a]) dim[0] *= legCosts[a];
    }
    else {
        dim = new double[2]; // *
        dim[0] = 1;
        dim[1] = 0;
        for (mwIndex a=0;a<numleglabels;a++) if (freelegs[a]) {
            dim[0] *= legCosts[a];
            dim[1] += legCosts[a+numleglabels];
        }
    }
    return dim;
}

inline double* dimSquared(const double *dim, const mwIndex costType)
{
    double *newdim = new double[mwIndex(costType)]; // *
    newdim[0] = dim[0]*dim[0];
    if (costType==2) newdim[1] = dim[1]*2;
    return newdim;
}

inline bool isgreaterthan_sd(const double *dim1, const double *dim2, const mwIndex costType)
{
    // Compares two single-index costs
    if (costType==1) return dim1[0]>dim2[0];
    if (dim1[1] > dim2[1]) return true;
    if (dim1[1] < dim2[1]) return false;
    if (dim1[0] > dim2[0]) return true;
    return false;
}

inline double* getbuildcost(const bool *freelegs, const bool *commonlegs, const double *legCosts, const mwIndex costType, const double oldmuCap, const double muCap, bool isnew,
        const double *cost1, const double *cost2, const mwIndex cost1len, const mwIndex cost2len, const mwIndex numleglabels, double &rtnmuCap, bool &isOK, mwIndex &newCostLen)
{
    // Get fusion cost
    bool *allLegs = new bool[numleglabels]; // OK
    for (mwIndex x=0;x<numleglabels;x++) allLegs[x] = freelegs[x] || commonlegs[x];
    if (costType==1) {
        double cost;
        cost = 1;
        for (mwIndex x=0;x<numleglabels;x++) if (allLegs[x]) cost *= legCosts[x];
        cost += cost1[0] + cost2[0];
        
        // ### Is cost too high (>muCap)?
        if (cost > muCap) {
            isOK = false; rtnmuCap = cost; delete[] allLegs; allLegs=NULL; return NULL;
        }
        // ### Is cost too low (not made from new objects, and <=oldmuCap: This construction has been done before)
        if (!isnew && cost <= oldmuCap) {
            isOK = false; rtnmuCap = 0; delete[] allLegs; allLegs=NULL; return NULL;
        } 
        // ### Done checking bounds on cost
        
        // Return new cost
        double *newCost = new double[1]; // *
        newCostLen = 1;
        newCost[0] = cost;
        delete[] allLegs; allLegs=NULL;
        return newCost;
    }
    else { // costType==2
        double fusionpower = 0;
        for (mwIndex x=0;x<numleglabels;x++) if (allLegs[x]) fusionpower += legCosts[numleglabels+x];
        
        // ### Is cost too high (>muCap)?
        if (fusionpower > muCap) {isOK = false; rtnmuCap = fusionpower; delete[] allLegs; allLegs=NULL; return NULL;}
        // ### Is cost too low (not made from new objects, and <=oldmuCap: This construction has been done before)
        if (!isnew && fusionpower <= oldmuCap) {isOK = false; rtnmuCap = 0; delete[] allLegs; allLegs=NULL; return NULL;}
        // ### Done checking bounds on cost
        
        // If cost OK, determine total cost of construction
        newCostLen = (cost1len > cost2len) ? cost1len : cost2len;
        if (newCostLen < fusionpower+1) newCostLen = fusionpower+1;
        double *newCost = new double[newCostLen]; // *
        for (mwIndex x=0;x<newCostLen;x++) {
            newCost[x] = 0;
            if (x<cost1len) newCost[x] += cost1[x];
            if (x<cost2len) newCost[x] += cost2[x];
        }
        double factor = 1;
        for (mwIndex x=0;x<numleglabels;x++) if (allLegs[x]) factor *= legCosts[x];
        newCost[(mwIndex)fusionpower] += factor;
        delete[] allLegs; allLegs=NULL;
        return newCost;
    }
}

inline bool islessthan(const double *cost1,const double *cost2,const mwIndex len1,const mwIndex len2,const mwIndex costType)
{
    // Compares two full network costs
    if (costType==1) return cost1[0] < cost2[0];
    if (len1 < len2) return true;
    if (len1 > len2) return false;
    for (mwIndex x=len2-1;x>=0;x--) {
        if (cost1[x] < cost2[x]) return true;
        if (cost1[x] > cost2[x]) return false;
    }
    return false;
}

inline void displaycostandsequence(const objectData *ptr, const double costtype, const mxArray* tracedindices, const mxArray *posindices)
{
    const double *cost = ptr->cost;
    const mwIndex costlen = ptr->costlen;
    const mwIndex *legsequence = ptr->sequence;
    const mwIndex legsequencelen = ptr->sequencelen;
    
    mexPrintf("\nSequence:      ");
    const int *tracedindices_re = (int*)mxGetPr(tracedindices);
    const double *posindices_re = mxGetPr(posindices);
    for (mwIndex x=0;x<mxGetNumberOfElements(tracedindices);x++) mexPrintf(" %g",posindices_re[tracedindices_re[x]-1]);
    for (mwIndex x=0;x<legsequencelen;x++) {
		if (legsequence[x]==0) mexPrintf(" 0");
		else {
            mwIndex legseq_t = legsequence[x];
            for (mwIndex y=0;y<mxGetNumberOfElements(tracedindices);y++) {
                if (tracedindices_re[y]<=legseq_t) legseq_t++;
                else break;
            }
            mexPrintf(" %g",posindices_re[legseq_t-1]);
        }
	}
    mexPrintf("\n");
    mexPrintf("Cost:           ");
    if (costtype==2) {
        for (mwIndex x=costlen-1;x>0;x--)
            mexPrintf("%gX^%d + ",cost[x],x);
        mexPrintf("%gX^0",cost[0]);
    }
    else mexPrintf("%g",cost[0]);
    if (mxGetNumberOfElements(tracedindices)!=0) mexPrintf(" + tracing costs");
    mexPrintf("\n");
    pause(0.01);
}

inline void addToTensorXlist(tensorXlist *&Xlistroot, const bool *freelegs, const mwIndex numleglabels, const double *newallin, const mwIndex costType, const bool isnew, const bool oldMayHaveEntry)
{
    // Constructed a new tensor for tensorXlist, or a known tensor at same or better cost.
    // If legs exactly match a known non-provisional entry, consider as a possible tighter bound on allIn.
    // Otherwise, consider for provisional list if not made redundant by any non-provisional entries.
    // Add to provisional list if not made redundant by any non-provisional entries.
    bool consider = true;
    tensorXlist *Xlistptr = Xlistroot;
    bool done = false;
    while (!done && Xlistptr!=NULL) {
        if (Xlistptr->flag==2) done = true;
        else {
            if (equalflags(Xlistptr->legs,freelegs,numleglabels)) {
                // If legs exactly match a non-provisional entry, update value of allIn.
                // If allIn for this entry just got increased, associated constraints have been relaxed. Flag this updated entry as provisional to trigger another pass.
                if (isgreaterthan_sd(newallin,Xlistptr->allIn,costType)) {
                    for (mwIndex x=0;x<costType;x++) Xlistptr->allIn[x] = newallin[x];
                    Xlistptr->flag = 2;
                    // Flags are always in ascending order: Move re-flagged entry to end of list
                    if (Xlistptr == Xlistroot) {
                        if (Xlistroot->next != NULL) {
                            Xlistroot = Xlistroot->next;
                            Xlistroot->prev = NULL;
                        }
                    }
                    else {
                        Xlistptr->prev->next = Xlistptr->next;
                        if (Xlistptr->next!=NULL) Xlistptr->next->prev = Xlistptr->prev;
                    }
                    if (Xlistptr != Xlistroot) {
                        tensorXlist *Xlistptr2 = Xlistroot;
                        while (Xlistptr2->next!=NULL) Xlistptr2 = Xlistptr2->next;
                        Xlistptr2->next = Xlistptr;
                        Xlistptr->prev = Xlistptr2;
                        Xlistptr->next = NULL;
                    }
                }
                else {
                    for (mwIndex x=0;x<costType;x++) Xlistptr->allIn[x] = newallin[x];
                }
                consider = false;
                break;
            }
            else {
                // Check to see if made redundant by existing non-provisional entry
                if (allAofB(Xlistptr->legs,freelegs,numleglabels) && !isgreaterthan_sd(newallin,Xlistptr->allIn,costType)) {
                    // All legs in freelegs are in overlap with tensorXlegs{a}, and dimension of absorbed tensor is not greater: Proposed new entry is redundant.
                    // (Greater dimension is allowed as it would mean more permissive bounds for a subset of legs)
                    consider = false;
                    if (!isnew && oldMayHaveEntry) {
                        // This tensor: Excluded from list.
                        // Previous, higer-cost contraction sequence may have successfully made an entry. Remove it.
                        removeFromTensorXlist(Xlistroot,freelegs,numleglabels);
                    }
                    break;
                }
            }
            Xlistptr = Xlistptr->next;
        }
    }
    if (consider) {
        // Tensor not excluded after comparison with non-provisional entries:
        // If legs exactly match another provisional entry, new submission is a cheaper sequence so keep only the new value of allIn.
        // (This is the same subnetwork but with a better cost, and possibly a different value of \xi_C.)
        if (!isnew && oldMayHaveEntry) {
            while (Xlistptr!=NULL) {
                if (Xlistptr->flag==2) {
                    if (equalflags(Xlistptr->legs,freelegs,numleglabels)) {
                        for (mwIndex x=0;x<costType;x++) Xlistptr->allIn[x] = newallin[x];
                        consider = false;
                        break;
                    }
                }
                Xlistptr = Xlistptr->next;
            }
        }
        // Otherwise: This is a new tensorXlist entry
        if (consider) {
            if (Xlistptr==NULL) Xlistptr = Xlistroot;
            while (Xlistptr->next!=NULL) Xlistptr = Xlistptr->next;
            Xlistptr->next = new tensorXlist; // OK
            Xlistptr->next->prev = Xlistptr;
            Xlistptr = Xlistptr->next;
            Xlistptr->legs = new bool[numleglabels]; // OK
            for (mwIndex x=0;x<numleglabels;x++) Xlistptr->legs[x] = freelegs[x];
            Xlistptr->allIn = new double[mwIndex(costType)]; // OK
            for (mwIndex x=0;x<costType;x++) Xlistptr->allIn[x] = newallin[x];
            // Xlistptr->flag defaults to 2
        }
    }
}

inline void updateTensorXlist(tensorXlist *&Xlistroot, const bool *freelegs, const mwIndex numleglabels, const double *newallin, const mwIndex costType)
{
    // Found new contraction sequence for an existing entry, but with a smaller value of allIn. Update the stored value to this value.
    // (This is the same subnetwork, but the \xi_C constraint has been tightened.)
    tensorXlist *Xlistptr = Xlistroot;
    while (Xlistptr!=NULL) {
        if (equalflags(Xlistptr->legs,freelegs,numleglabels)) {
            for (mwIndex x=0;x<costType;x++) Xlistptr->allIn[x] = newallin[x];
            break;
        }
        Xlistptr = Xlistptr->next;
    }
    // If the original entry was provisional and redundant, this one is too. No match will be found (as the redundant tensor was never recorded), and that's OK.
    // If the original entry was not provisional and redundant, it has now been updated.
}

inline void removeFromTensorXlist(tensorXlist *&Xlistroot, const bool *freelegs, const mwIndex numleglabels)
{
    // Remove this tensor from tensorXlist.
    // Cannot restrict checking just to provisional portion, because might need to remove a confirmed tensor's data from the list in the scenario described in footnote 2.
    tensorXlist *Xlistptr = Xlistroot;
    while (Xlistptr!=NULL) {
        if (equalflags(Xlistptr->legs,freelegs,numleglabels)) {
            if (Xlistptr==Xlistroot) {
                Xlistroot = Xlistroot->next;
                Xlistroot->prev = NULL;
                delete Xlistptr;
            }
            else {
                Xlistptr->prev->next = Xlistptr->next;
                if (Xlistptr->next!=NULL) Xlistptr->next->prev = Xlistptr->prev;
                delete Xlistptr;
            }
            break;
        }
        Xlistptr = Xlistptr->next;
    }
}

inline void mergeTensorXlist(tensorXlist *&Xlistroot, const mwIndex costType, const mwIndex numleglabels)
{
    // Merge provisional tensors into main list and decrease all nonzero flags by one.
    // For each provisional entry in turn, check if it makes redundant or is made redundant by any other provisional entries, and if it makes redundant any
    // non-provisional entries. If so, delete the redundant entries.
    
    tensorXlist *Xlistptr = Xlistroot;
    tensorXlist *firstprovtensor = NULL;

    // Decrease non-zero tensorXflags and find first provisional tensor
    while (Xlistptr!=NULL) {
        if (Xlistptr->flag==2 && firstprovtensor==NULL) firstprovtensor = Xlistptr;
        if (Xlistptr->flag > 0) Xlistptr->flag--;
        Xlistptr = Xlistptr->next;
    }
    
    // Permit provisional and non-provisional tensors to be eliminated by other provisional tensors (elimination of provisional tensors by non-provisional tensors was done earlier)
    tensorXlist *ptrA = Xlistroot;
    while (ptrA!=NULL) {
        bool advance = true;
        tensorXlist *ptrB = firstprovtensor;
        while (ptrB!=NULL) {
            if (ptrA!=ptrB) {
                if (allAofB(ptrB->legs,ptrA->legs,numleglabels) && !isgreaterthan_sd(ptrA->allIn,ptrB->allIn,costType)) { // All legs in tensorXlegs{a} overlap with tensorXlegs{b}, and dimension of absorbed tensor is not greater: tensorXlegs{a} is redundant.
                    if (ptrA==Xlistroot) {
                        Xlistroot = Xlistroot->next;
                        Xlistroot->prev = NULL;
                        delete ptrA;
                        ptrA = Xlistroot;
                        advance = false;
                    }
                    else {
                        ptrA->prev->next = ptrA->next;
                        if (ptrA->next!=NULL) ptrA->next->prev = ptrA->prev;
                        tensorXlist *t = ptrA->next;
                        delete ptrA;
                        ptrA = t;
                        advance = false;
                    }
                    break;
                }
            }
            ptrB = ptrB->next;
        }
        if (advance) ptrA = ptrA->next;
    }
}

inline bool allundertwo(tensorXlist *list)
{
    while (list!=NULL) {
        if (list->flag==2) return false;
        list = list->next;
    }
    return true;
}
