
/*
 *
 * TODO: Fix the following bug: When generating the list of all normalized vertices, it's not enough to cycle through all partitions, we also have to ensure that the normalized vertices are maximal with respect to all permutations of the skipped vertices.
 *       This means we have to cycle through all partitions, generate the normpts, use npts2pts to get the points, cycle through all permutations, generate the normpts from the permuted points and then compare the OG normpts with the perm normpts to see if they are maximal.
 * TODO: there seems to be a bug in the generation of all normed verts, even w/o maximizing. We don't get all combos.
 * TODO: there seems to be a bug in the computation of the normalized form. For some verts we compute a form which is not maximal.
 * 
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

typedef unsigned long ulong;
typedef unsigned char uchar;
typedef struct vertstruct vert;
typedef struct cubestruct cube;
typedef struct loosepathstruct loosepath;
typedef struct loosepathslenstruct loosepathslen;
typedef struct loosepathsstruct loosepaths;
typedef struct normptsstruct normpts;
typedef struct normptscntstruct normptscnt;
typedef struct permsstruct perms;

struct vertstruct {
    uchar *label;
    uchar covered;
};

struct cubestruct {
    uchar k;
    uchar d;
    vert *v;
    ulong vcnt;
    ulong ecnt;
    ulong lhplen;
};

struct normptsstruct {
    uchar val;
    uchar max;
    uchar cntpos;
    uchar chcnt;
    normpts **children;
};

struct normptscntstruct {
    normpts *npts;
    uchar **nverts;
    ulong *perm;
    uchar ismax;
    normptscnt *max;
    ulong maxid;
    ulong cnt;
};

struct loosepathstruct {
    ulong len;
    vert **v;
    vert **skipped;
    normptscnt *npc;
    loosepath *next;
};

struct loosepathslenstruct {
    ulong cnt,lene,lenv,lenskipped,nptsleafcnt,nptscnt,nptsmaxcnt,nptscovered;
    normptscnt **npts;
    loosepath **lps;
};

struct loosepathsstruct {
    ulong cnt;
    loosepath *first;
    loosepath *last;
    ulong lmin;
    ulong lmax;
    ulong ldel;
    loosepathslen *lpl;
};

struct permsstruct {
    ulong cntels;
    ulong cntps;
    ulong **ps;
};

ulong powuluc(ulong a, uchar b) {
    ulong res;
    uchar i;
    res=1;
    for(i=0;i<b;i++) {
        res=res*a;
    }
    return res;
}

ulong powulul(ulong a, ulong b) {
    ulong res;
    ulong i;
    res=1;
    for(i=0;i<b;i++) {
        res=res*a;
    }
    return res;
}

ulong factorial(ulong a) {
    ulong res;
    ulong b;
    res=1;
    for(b=2;b<=a;b++) {
        res*=b;
    }
    return res;
}

ulong binomial(ulong a,ulong b) {
    double x;
    ulong i,a2,b2;
    if(b>a) return 0;
    if(b>b-a) b=b-a;
    x=1; a2=a; b2=b;
    for(i=0;i<b;i++) {
        x=x*a2/b2;
        a2--;b2--;
    }
    return (ulong) rint(x);
}

ulong pt2ind(uchar len, uchar *a, uchar b) {
    ulong res, p;
    uchar i;
    p=1;
    res=0;
    for(i=0;i<len;i++) {
        res+=a[i]*p;
        p*=b;
    }
    return res;
}

ulong pt2indul(ulong len, ulong *a, ulong b) {
    ulong res, p;
    ulong i;
    p=1;
    res=0;
    for(i=0;i<len;i++) {
        res+=a[i]*p;
        p*=b;
    }
    return res;
}

void ind2pt(ulong ind, uchar b, uchar len, uchar *a) {
    ulong x,y;
    uchar i;
    x=ind;
    for(i=0;i<len;i++) {
        y=x/b;
        a[i]=x-b*y;
        x=y;
    }
}

void ind2ptul(ulong ind, ulong b, ulong len, ulong *a) {
    ulong x,y;
    ulong i;
    x=ind;
    for(i=0;i<len;i++) {
        y=x/b;
        a[i]=x-b*y;
        x=y;
    }
}

uchar isinj(ulong *a, ulong len) {
    ulong i,i2;
    for(i=0;i<len;i++) {
        for(i2=0;i2<i;i2++) {
            if(a[i]==a[i2]) return 0;
        }
    }
    return 1;
}

void initcube(cube *c) {
    ulong i;
    c->v=malloc(c->vcnt*sizeof(vert));
    for(i=0;i<c->vcnt;i++) {
        c->v[i].covered=0;
        c->v[i].label=malloc(c->d*sizeof(uchar));
        ind2pt(i,c->k,c->d,c->v[i].label);
    }
}

void freecube(cube *c) {
    free(c->v);
}

vert *label2vert(cube *c,uchar *label) {    
    return c->v+pt2ind(c->d,label,c->k);
}

// TODO: find an efficient algorithm which initializes perms to embed symmetric groups S1 setle ... setle S{cntels}, meaning
//       ps[0]=identity (-> S1), ps[1]=21 (rest id), ps[2]=312, ps[3]=321, ps[4]=132, ps[5]=231, ps[6]=4123, ...
perms *initperms(ulong cntels) {
    ulong cntpts,i,i2;
    uchar good;
    perms *ps=malloc(sizeof(perms));
    ps->cntels=cntels;
    ps->cntps=factorial(cntels);
    ps->ps=malloc(ps->cntps*sizeof(ulong *));
    if(cntels==0) {
        ps->ps[0]=NULL;
        return ps;
    }
    for(i=0;i<ps->cntps;i++) ps->ps[i]=malloc(cntels*sizeof(ulong));
    cntpts=powulul(cntels,cntels);
    i=0;
    for(i2=0;i<cntpts;i2++) {
        ind2ptul(i2,cntels,cntels,ps->ps[i]);
        if(isinj(ps->ps[i],cntels)) {
            if(i==ps->cntps-1) return ps;
            i+=1;
        }
    }
    // we never end up here b/c the last pt is all-m, which is not a perm
    return ps;
}

void freeperms(perms *ps) {
    ulong i;
    if(ps->cntels>0) {
        for(i=0;i<ps->cntps;i++) free(ps->ps[i]);
    }
    free(ps->ps);
    free(ps);
}

normpts *initnormpt(cube *c,uchar val, uchar max, ulong h) {
    uchar i;
    normpts *npts=malloc(sizeof(normpts));
    npts->cntpos=0;
    npts->val=val;
    npts->max=(val>max)?val:max;
    if(h<2) {
        npts->chcnt=0;
        npts->children=NULL;
        return npts;
    }
    // max < k-1 => we get all values from 0 to max+1 for children, so chcnt=max+2
    npts->chcnt=(npts->max+2>c->k)?c->k:npts->max+2;
    npts->children=malloc(npts->chcnt*sizeof(normpts *));
    for(i=0;i<npts->chcnt;i++) {
        npts->children[i]=initnormpt(c,i,npts->max,h-1);
    }
    return npts;
}

void freenormpt(normpts *npts) {
    uchar m;
    for(m=0;m<npts->chcnt;m++) {
        freenormpt(npts->children[m]);
    }
    if(npts->chcnt>0) free(npts->children);
    free(npts);
}

void nptsreset(normpts *npts) {
    uchar i;
    for(i=0;i<npts->chcnt;i++) nptsreset(npts->children[i]);
    npts->cntpos=0;
}

char nptscompare(normpts *lhs,normpts *rhs) {
    uchar i;
    char cmp;
    if(lhs->chcnt==0) {
        if(lhs->cntpos>rhs->cntpos) return 1;
        if(lhs->cntpos<rhs->cntpos) return -1;
        return 0;
    }
    for(i=0;i<lhs->chcnt;i++) {
        cmp=nptscompare(lhs->children[i],rhs->children[i]);
        if(cmp!=0) return cmp;
    }
    return 0;
}

void nptscopy(normpts *src,normpts *tgt) {
    uchar i;
    tgt->cntpos=src->cntpos;
    for(i=0;i<src->chcnt;i++) nptscopy(src->children[i],tgt->children[i]);
}

ulong npts2pts(normpts *npts, ulong h, uchar **pts,ulong ind,ulong ccnt) {
    uchar i;
    ulong j;
    if(ind>=ccnt) return ind;
    pts[h][ind]=npts->val;
    if(npts->chcnt==0 && npts->cntpos>0) {
        for(i=1;i<npts->cntpos;i++) {
            for(j=0;j<=h;j++) pts[j][ind+i]=pts[j][ind];
        }
        ind+=npts->cntpos;
        if(ind<ccnt) for(j=0;j<=h;j++) pts[j][ind]=pts[j][ind-1];
        return ind;
    }
    for(i=0;i<npts->chcnt;i++) {
        ind=npts2pts(npts->children[i],h+1,pts,ind,ccnt);
    }
    return ind;
}

ulong nptsleafcnt(normpts *npts) {
    ulong res,i;
    if(npts->chcnt==0) return 1;
    res=0;
    for(i=0;i<npts->chcnt;i++) {
        res+=nptsleafcnt(npts->children[i]);
    }
    return res;
}

ulong nptspopulateRec(normpts *npts, uchar *vals, ulong ind) {
    uchar i;
    if(npts->chcnt==0) {
        npts->cntpos=vals[ind];
        return ind+1;
    }
    npts->cntpos=0;
    for(i=0;i<npts->chcnt;i++) {
        ind=nptspopulateRec(npts->children[i],vals,ind);
        npts->cntpos+=npts->children[i]->cntpos;
    }
    return ind;
}

void nptspopulate(normpts *npts, uchar *vals) {
    nptspopulateRec(npts,vals,0);
}

void normpts2pts(normpts *npts, uchar **pts) {
    npts2pts(npts,0,pts,0,npts->cntpos);
}

void pts2normpts(cube *c, loosepathslen *lpl, uchar **pts, normpts *npts) {
    uchar max,coord,oldval;
    ulong pt,pt2;
    normpts *ptr;
    nptsreset(npts);
    for(coord=0;coord<c->d;coord++) {
        max=0;
        for(pt=0;pt<lpl->lenskipped+2;pt++) {
            if(pts[pt][coord]>max) max++;
            if(pts[pt][coord]>max) {
                oldval=pts[pt][coord];
                pts[pt][coord]=max;
                for(pt2=pt+1;pt2<lpl->lenskipped+2;pt2++) {
                    if(pts[pt2][coord]==max) pts[pt2][coord]=oldval;
                    else if(pts[pt2][coord]==oldval) pts[pt2][coord]=max;
                }
            }
        }
    }
    for(coord=0;coord<c->d;coord++) {
        ptr=npts;
        ptr->cntpos++;
        for(pt=1;pt<lpl->lenskipped+2;pt++) {
            ptr=ptr->children[pts[pt][coord]];
            ptr->cntpos++;
        }
    }
}

normpts *initnormpts(cube *c, ulong scnt) {
    return initnormpt(c,0,0,scnt+2);
}

uchar pointsdistinct(uchar ccnt,ulong ptcnt,uchar **pts) {
    ulong i,i2;
    uchar j,equal;
    for(i=0;i<ptcnt;i++) {
        for(i2=0;i2<i;i2++) {
            equal=1;
            for(j=0;j<ccnt;j++) {
                if(pts[i][j]!=pts[i2][j]) equal=0;
            }
            if(equal) return 0;
        }
    }
    return 1;
}

normptscnt *initnptscnt(cube *c,loosepathslen *lpl) {
    normptscnt *npc;
    ulong i;
    npc=malloc(sizeof(normptscnt));
    npc->npts=initnormpts(c,lpl->lenskipped);
    npc->nverts=malloc((lpl->lenskipped+2)*sizeof(uchar *));
    for(i=0;i<lpl->lenskipped+2;i++) npc->nverts[i]=malloc(c->d*sizeof(uchar));
    npc->cnt=0;
    npc->ismax=0;
    npc->max=NULL;
    npc->perm=malloc(lpl->lenskipped*sizeof(long));
    npc->maxid=0;
    return npc;
}

void freenptscnt(loosepathslen *lpl,normptscnt *npc) {
    ulong i;
    free(npc->perm);
    for(i=0;i<lpl->lenskipped+2;i++) free(npc->nverts[i]);
    free(npc->nverts);
    freenormpt(npc->npts);
    free(npc);
}

ulong filllplnpts(cube *c, loosepathslen *lpl, uchar *vals, uchar coord,uchar sum,ulong ind) {
    uchar i;
    if(coord==lpl->nptsleafcnt-1) {
        vals[coord]=sum;
        nptspopulate(lpl->npts[ind]->npts,vals);
        normpts2pts(lpl->npts[ind]->npts,lpl->npts[ind]->nverts);
        if(!pointsdistinct(c->d,lpl->lenskipped+2,lpl->npts[ind]->nverts)) return ind;
        ind=ind+1;
        lpl->npts[ind]=initnptscnt(c,lpl);
        return ind;
    }
    for(i=0;i<=sum;i++) {
        vals[coord]=i;
        ind=filllplnpts(c,lpl,vals,coord+1,sum-i,ind);
    }
    return ind;
}

void normptscntfindmax(cube *c, loosepathslen *lpl) {
    perms *ps;
    ulong npc,npc2,p,pt;
    uchar **nvals, coord;
    normpts *npts, *npts2;
    lpl->nptsmaxcnt=0;
    ps=initperms(lpl->lenskipped);
    nvals=malloc((lpl->lenskipped+2)*sizeof(uchar *));
    for(pt=0;pt<lpl->lenskipped+2;pt++) nvals[pt]=malloc(c->d*sizeof(uchar));
    for(coord=0;coord<c->d;coord++) nvals[0][coord]=0;
    npts=initnormpts(c,lpl->lenskipped);
    npts2=initnormpts(c,lpl->lenskipped);
    for(npc=0;npc<lpl->nptscnt;npc++) {
        for(coord=0;coord<c->d;coord++) nvals[1][coord]=lpl->npts[npc]->nverts[1][coord];
        for(p=0;p<ps->cntps;p++) {
            for(pt=0;pt<lpl->lenskipped;pt++) {
                for(coord=0;coord<c->d;coord++) nvals[pt+2][coord]=lpl->npts[npc]->nverts[ps->ps[p][pt]+2][coord];
            }
            pts2normpts(c,lpl,nvals,npts);
            if(nptscompare(npts,npts2)>0) {
                nptscopy(npts,npts2);
                for(pt=0;pt<lpl->lenskipped;pt++) lpl->npts[npc]->perm[pt]=ps->ps[p][pt];
            }
        }
        npc2=0;
        while(nptscompare(npts2,lpl->npts[npc2]->npts)!=0) npc2++;
        lpl->npts[npc]->max=lpl->npts[npc2];
        if(npc2==npc) {
            lpl->npts[npc]->ismax=1;
            lpl->nptsmaxcnt++;
        }
        nptsreset(npts2);
    }
    freeperms(ps);
    for(pt=0;pt<lpl->lenskipped+2;pt++) free(nvals[pt]);
    free(nvals);
    freenormpt(npts);
    freenormpt(npts2);
    npc2=0;
    for(npc=0;npc<lpl->nptscnt;npc++) {
        if(lpl->npts[npc]->ismax) {
            lpl->npts[npc]->maxid=npc2;
            npc2++;
        }
    }
}

void initlplnpts(cube *c,loosepathslen *lpl) {
    ulong i;
    ulong len;
    uchar *vals;
    vals=malloc(lpl->nptsleafcnt*sizeof(uchar));
    len=binomial(c->d+lpl->nptsleafcnt-1,c->d);
    lpl->npts=malloc(len*sizeof(normptscnt *));
    lpl->nptscnt=0;
    lpl->npts[0]=initnptscnt(c,lpl);
    lpl->nptscnt=filllplnpts(c,lpl,vals,0,c->d,0);
    freenptscnt(lpl,lpl->npts[lpl->nptscnt]);
    free(vals);
    normptscntfindmax(c,lpl);
}

void initlps(cube *c, loosepaths *lps) {
    ulong i;
    normpts *npts;
    lps->first=malloc(sizeof(loosepath));
    lps->last=lps->first;
    lps->cnt=0;
    lps->ldel=lps->lmax-lps->lmin+1;
    lps->lpl=malloc(lps->ldel*sizeof(loosepathslen));
    for(i=0;i<lps->ldel;i++) {
        lps->lpl[i].cnt=0;
        lps->lpl[i].lene=lps->lmin+i;
        lps->lpl[i].lenv=1+lps->lpl[i].lene*(c->k-1);
        lps->lpl[i].lenskipped=(c->lhplen-lps->lpl[i].lene)*(c->k-1);
        npts=initnormpts(c,lps->lpl[i].lenskipped);
        lps->lpl[i].nptsleafcnt=nptsleafcnt(npts);
        freenormpt(npts);
        initlplnpts(c,lps->lpl+i);
        lps->lpl[i].nptscovered=0;
    }
}

void freelps(loosepaths *lps) {
    ulong i,j;
    loosepath *lp, *lp2;
    lp=lps->first;
    lp2=lp->next;
    free(lp);
    for(i=0;i<lps->cnt;i++) {
        lp=lp2;
        lp2=lp->next;
        free(lp->v);
        free(lp->skipped);
        free(lp);
    }
    for(i=0;i<lps->ldel;i++) {
        free(lps->lpl[i].lps);
        for(j=0;j<lps->lpl[i].nptscnt;j++) freenptscnt(lps->lpl+i,lps->lpl[i].npts[j]);
        free(lps->lpl[i].npts);
    }
    free(lps->lpl);
}

void normalizelp(cube *c, loosepathslen *lpl, ulong i, normpts *npts, uchar **nvals) {
    loosepath *lp;
    uchar coord;
    ulong pt, npc;
    lp=lpl->lps[i];
    for(coord=0;coord<c->d;coord++) nvals[1][coord]=lp->v[lpl->lenv-1]->label[coord];
    for(pt=0;pt<lpl->lenskipped;pt++) {
        for(coord=0;coord<c->d;coord++) nvals[2+pt][coord]=lp->skipped[pt]->label[coord];
    }
    pts2normpts(c,lpl,nvals,npts);
}

void normalizelpl(cube *c, loosepathslen *lpl) {
    ulong i,j,is;
    loosepath *lp;
    normpts *npts=initnormpts(c,lpl->lenskipped);
    uchar **nvals=malloc((lpl->lenskipped+2)*sizeof(ulong *));
    for(i=0;i<lpl->lenskipped+2;i++) nvals[i]=malloc(c->d*sizeof(uchar));
    for(i=0;i<c->d;i++) nvals[0][i]=0;
    for(i=0;i<lpl->cnt;i++) {
        lp=lpl->lps[i];
        normalizelp(c,lpl,i,npts,nvals);
        j=0;
        while(nptscompare(lpl->npts[j]->npts,npts)!=0) j++;
        lp->npc=lpl->npts[j];
        lpl->npts[j]->max->cnt++;
    }
    for(i=0;i<lpl->lenskipped+2;i++) free(nvals[i]);
    free(nvals);
    freenormpt(npts);
    for(i=0;i<lpl->nptscnt;i++) {
        if(lpl->npts[i]->cnt>0) lpl->nptscovered++;
    }
}

/*
 *
 * EXPOSITION
 * All coord & val per coord perms break the normalization of LHPs in lps
 * Reason: for given LHP and nontriv perms pi the result LHP' must violate one of the minmal 1st-exploration constraints
 *    For coord perm exists c with pic(c)<c, then in LHP' pic(c) is still the c-th coord to be explored, but should be the pic(c)-th coord to be explored
 *    For vals it's analogous, take c & x with piv{c}(x)<x, then piv{c}(x) is still only the x-th val to be explored in LHP'
 * So, no two LHPs in lps are equivalent & no LHP is invariant under coord & val per coord perms
 * Other LPs in lps can be invariant under coord & val per coord perms (think of the 3 unique len <3 LPs)
 * However, two LPs still cannot be equivalent
 * Reason: If you apply perms to an LHP such that the result LHP' differs from LHP, then LHP' necessarily violates normalization (as before)
 * So, if we apply perms to the LPs in lps, we maintain uniqueness of the LPs
 * 
 * EXPLANATION
 * normalizelps computes normalizations of start w0, end w1 & skipped vertices s0,...,s{m-1}
 * For the normalization we can use the following types of permutations
 * 1. permutations of coordinates (allowed b/c this gives equivalent LP)
 * 2. permutations of values for given coordinate keeping 0 fixed (allowed b/c this gives equivalent LP, enforcing that w0 stays the origin)
 * 3. permutations of the skipped verts, i.e. perms pi:{0,...,m-1}->{0,...,m-1} (allowed b/c skipped verts model a set {s0,...,s{m-1})
 * For given w0, w1, s0,...,s{m-1} we can use type 1 & 2 perms to construct matrices in {0,...,k-1}^{d x m+2} of the form
 * 0 0 0 0
 * 0 0 0 0
 * 0 0 0 1
 * 0 0 0 1
 * 0 0 1 0
 * 0 0 1 0
 * 0 0 1 1
 * 0 0 1 1
 * 0 0 1 2
 * 0 0 1 2
 * 0 1 0 0
 * 0 1 0 0
 * 0 1 0 1
 * 0 1 0 1
 * 0 1 0 2
 * 0 1 0 2
 * 0 1 1 0
 * 0 1 1 0
 * 0 1 1 1
 * 0 1 1 1
 * 0 1 1 2
 * 0 1 1 2
 * 0 1 2 0
 * 0 1 2 0
 * 0 1 2 1
 * 0 1 2 1
 * 0 1 2 2
 * 0 1 2 2
 * 0 1 2 3
 * 0 1 2 3
 * A more efficient way to encode these matrices is a tree, where each node comes with a value, a max and a cnt
 * The root has val 0, max 0, cnt d; node with max has min(k,max+1) children; the x-th child has val x & max max(max,x)
 * For each node the cnt is the number of coords in the matrix such that the vals from left to right match the vals on the path from the root to the node.
 * So for the root cnt=d (only 0s), for node given by 00 cnt=# 0s in 1st col, for 01 cnt=# 1s in 1st col, for 000 cnt=#0s in 2nd col where 1st col=0 
 * We define an order on these trees using dfs with the counts (for a fixed structure, ie trees only differ in cnt's)
 * Using type 3 perms we can choose the maximal tree
*/
void normalizelps(cube *c, loosepaths *lps) {
    ulong i,j;
    for(i=0;i<lps->ldel;i++) {
        normalizelpl(c,lps->lpl+i);
    }
}

void printdata(cube *c, loosepaths *lps, char *fin) {
    char *fname;
    ulong x,i,i2,i3,len;
    loosepathslen *lpl;
    loosepath *lp;
    FILE *f;
    uchar **nverts;
    len=strlen(fin)-4;
    fname=malloc((len+11)*sizeof(char));
    strncpy(fname,fin,len);
    sprintf(fname+len,"_FINAL.csv");
    f=fopen(fname,"w");
    while(f==NULL) {
        printf("Could not open file %s. Press ENTER to try again or CTRL+C to abort.",fname);
        getchar();
        f=fopen(fname,"w");
    }
    printf("\t\t\t\t\r");
    fprintf(f,"LOOSE PATHS IN THE CUBE HYPERGRAPH\n\n\nEXPOSITION AND DEFINITIONS\n\n");
    fprintf(f,"Cube Hypergraph,\"k-uniform d-regular hypergraph H(k,d) with vertices X^d, X={0,...,k-1}, and hyperedges {v_x:x in X}, where (v_x)_c=x, u=v_x o f_c and f_c:[d-1]->[d]\\{c} is the enumeration, for all u in [k]^{d-1} and c in [d]\"\n");
    fprintf(f,"Loose Hamilton path length (=hyperedge count),\"m=(k^d-1)/(k-1)\"\n");
    fprintf(f,"Minimum loose path (LP) length,\"m_min=m-dmin\"\n");
    fprintf(f,"Maximum LP length,\"m_max=m-dmax\"\n");
    fprintf(f,"Ordered LP,\"A bijection v:[k^d]->[k]^d such that e_i={v_{1+(i-1)(k-1)},...,v_{1+i(k-1)}} is a hyperedge in H(k,d) for all i in {1,...,m'} with m' in [m_min,m_max]\"\n,\"The starting point is v_1 and the end point is v_{1+m'(k-1)}\"\n");
    fprintf(f,"Directed LP,\"Directed LP for ordered LP v is (u_0,e_1,u_1,e_2,...,u_{m'-1},e_m',u_m'), where u_0=v_1, e_i from above and u_i=v_{1+i(k-1)}\"\n,\"Starting point is u_0 and end point is u_m'\"\n");
    fprintf(f,"Undirected LP,\"Undirected LP for directed LP (u_0,e_1,...,e_m',u_m') is {(u_0,e_1,...,e_m',u_m'),(u_m',e_m',...,e_1,u_0)}\"\n,\"End points are u_0 and u_m'\"\n");
    fprintf(f,"Graph isomorphisms,\"For each vertex v in H(k,d) fix an automorphism f_v:H(k,d)->H(k,d) that takes v to the origin\"\n,\"Any directed LP with starting point v is taken to a directed LP with the origin as starting point under f_v\"\n");
    fprintf(f,"Type 1 LPs,\"Directed LP starting at the origin with the following properties\"\n,\"Let i_j be the step where the j-th 0 was changed for the first time, e.g. i_1=1 because we change a coordinate of the origin in the first step and i_2=2 because we have to change another coordinate in the second step. Let c_j be the unique coordinate in which u_{i_j} and u_{i_j-1} differ. Then we have c_j=j for all j in [d]\"\n,\"Let i_{c,j} be the step where u_i takes a value (in X) on coordinate c that has not appeared before (on coordinate c) for the j-th time. Let x_{c,j} be this new value. Then we have x_{c,j}=j for all j in [k-1] and all c in [d]\"\n");
    fprintf(f,"Type 2 LPs,\"A type 2 LP is a directed LP with the origin as starting point\"\n,\"Equivalently, a type 2 LP is given by a type 1 LP, a bijection [d]->[d] to fix the order in which the coordinates are explored and d bijections [k-1]->[k-1] to fix the order in which the values are explored, per coordinate\"\n");
    fprintf(f,"Type 3 LPs,\"directed LPs or equivalently given by a pair of a type 2 LP and a starting point v (using the fixed isomorphism f_v from above)\"\n");
    fprintf(f,"Skipped vertices,\"For a LP of any type, the skipped vertices S_0,...,S_{(m-m')(k-1)-1} are all vertices that do not appear in the LP\"");
    fprintf(f,"Normalized vertices,\"For a type 1 LP, we obtain representatives (w,s_0,...,s_{(m-m')(k-1)-1}) of (u_m',S_0,...,S_{(m-m')(k-1)-1}) under bijections [d]->[d] and [k-1]->[k-1]\"\n,\"We consider a canonical order on all representatives over all bijections and all LPs to obtain normalized vertices w,s_0,...,s_{(m-m')(k-1)-1}\"");
    fprintf(f,"\nRESULTS\n\nk,%u\nd,%u\n#vertices,%lu\n#edges,%lu\nm,%lu\nm_min,%lu\nm_max,%lu\n",c->k,c->d,c->vcnt,c->ecnt,c->lhplen,lps->lmin,lps->lmax);
    x=factorial(c->d)*powuluc(factorial(c->k-1),c->d);
    fprintf(f,"Type 2 LPs per Type 1 LP,%lu\n",x);
    fprintf(f,"Type 3 LPs per Type 2 LP,%lu\n",c->vcnt);
    x=x*c->vcnt;
    fprintf(f,"Type 3 LPs per Type 1 LP,%lu\n\n",x);
    fprintf(f,"TYPE 1 LPS (displayed as ordered LPs)\n\n");
    nverts=malloc((lps->lpl[0].lenskipped+2)*sizeof(uchar *));
    for(i=0;i<lps->lpl[0].lenskipped+2;i++) nverts[i]=malloc(c->d*sizeof(uchar));
    for(len=0;len<lps->ldel;len++) {
        lpl=lps->lpl+len;
        fprintf(f,"Edge count,%lu\nVertex Count,%lu\nSkipped Vertex Count,%lu\nNumber of LPs,%lu\nTotal Number of normalized vertices,%lu\nNumber of covered normalized vertices,%lu\n\n",lpl->lene,lpl->lenv,lpl->lenskipped,lpl->cnt,lpl->nptsmaxcnt,lpl->nptscovered);
        fprintf(f,"ID,w");
        for(i=0;i<lpl->lenskipped;i++) fprintf(f,",S[%lu]",i);
        fprintf(f,",LP count\n");
        for(i=0;i<lpl->nptscnt;i++) {
            if(lpl->npts[i]->ismax) {
                fprintf(f,"%lu,",lpl->npts[i]->maxid);
                for(i2=0;i2<c->d;i2++) fprintf(f,"%hhu",lpl->npts[i]->nverts[1][i2]);
                for(i2=0;i2<lpl->lenskipped;i2++) {
                    fprintf(f,",");
                    for(i3=0;i3<c->d;i3++) fprintf(f,"%hhu",lpl->npts[i]->nverts[2+i2][i3]);
                }
                fprintf(f,",%lu\n",lpl->npts[i]->cnt);
            }
        }
        fprintf(f,"\n");
        fprintf(f,"v[0]");
        for(i=1;i<lpl->lenv;i++) fprintf(f,",v[%lu]",i);
        for(i=0;i<lpl->lenskipped;i++) fprintf(f,",S[%lu]",i);
        fprintf(f,",nvind\n");
        for(i=0;i<lpl->cnt;i++) {
            lp=lpl->lps[i];
            for(i3=0;i3<c->d;i3++) fprintf(f,"%hhu",lp->v[0]->label[i3]);
            for(i2=1;i2<lpl->lenv;i2++) {
                fprintf(f,",");
                for(i3=0;i3<c->d;i3++) fprintf(f,"%hhu",lp->v[i2]->label[i3]);
            }
            for(i2=0;i2<lpl->lenskipped;i2++) {
                fprintf(f,",");
                for(i3=0;i3<c->d;i3++) fprintf(f,"%hhu",lp->skipped[lp->npc->perm[i2]]->label[i3]);
            }
            fprintf(f,",%lu\n",lp->npc->max->maxid);
        }
        fprintf(f,"\n");
    }
    fprintf(f,"\n");
    for(i=0;i<lps->lpl[0].lenskipped+1;i++) free(nverts[i]);
    free(nverts);
}

void process(char *fname) {
    cube c;
    loosepaths lps;
    loosepath *lp;
    FILE *f;
    long fsize;
    char *buf, *strcur;
    ulong i,j,ldel;
    uchar *label;
    f=fopen(fname,"r");
    fseek(f,0,SEEK_END);
    fsize=ftell(f);
    fseek(f,0,SEEK_SET);
    buf=malloc(fsize+1);
    fread(buf,1,fsize,f);
    fclose(f);
    strcur=strstr(buf,"\nk,")+1;
    sscanf(strcur,"k,%hhu\nd,%hhu\n#vertices,%lu\n#edges,%lu\nm,%lu\nm_min,%lu\nm_max,%lu\n",&(c.k),&(c.d),&(c.vcnt),&(c.ecnt),&(c.lhplen),&(lps.lmin),&(lps.lmax));
    initcube(&c);
    initlps(&c,&lps);
    strcur=strstr(strcur,"\nLength,")+1;
    strcur=strstr(strcur,"\n")+1;
    label=malloc(c.d*sizeof(uchar));
    // file complete: until empty line; file incomplete: as long as line completed
    while(strcur[0]>=48 && strstr(strcur,"\n")!=NULL) {
        lp=malloc(sizeof(loosepath));
        sscanf(strcur,"%lu",&(lp->len));
        strcur=strstr(strcur,",")+1;
        ldel=lp->len-lps.lmin;
        lp->v=malloc(lps.lpl[ldel].lenv*sizeof(vert *));
        for(i=0;i<lps.lpl[ldel].lenv;i++) {
            for(j=0;j<c.d;j++) label[j]=strcur[j]-48;
            lp->v[i]=label2vert(&c,label);
            strcur=strcur+c.d+1;
        }
        lp->next=NULL;
        lp->skipped=malloc(lps.lpl[ldel].lenskipped*sizeof(vert *));
        for(i=0;i<lps.lpl[ldel].lenv;i++) lp->v[i]->covered=1;
        j=0;
        for(i=0;i<c.vcnt;i++) {
            if(!c.v[i].covered) {
                lp->skipped[j]=c.v+i;
                j+=1;
            }
            c.v[i].covered=0;
        }
        lps.last->next=lp;
        lps.last=lp;
        lps.cnt+=1;
        lps.lpl[ldel].cnt+=1;
    }
    for(i=0;i<lps.ldel;i++) {
        lps.lpl[i].lps=malloc(lps.lpl[i].cnt*sizeof(loosepath *));
        lps.lpl[i].cnt=0;
    }
    lp=lps.first;
    for(i=0;i<lps.cnt;i++) {
        lp=lp->next;
        ldel=lp->len-lps.lmin;
        lps.lpl[ldel].lps[lps.lpl[ldel].cnt]=lp;
        lps.lpl[ldel].cnt+=1;
    }
    normalizelps(&c,&lps);
    printdata(&c,&lps,fname);
    free(label);
    freelps(&lps);
    freecube(&c);
    free(buf);
}

int main(int argc, char *argv[]) {
    uchar debug=0;
    char *fname, dfname[100];
    if(debug) {
        sprintf(dfname,"lhps_k3_d4_ALL.csv");
        fname=dfname;
    } else if(argc<2) {
        printf("Usage: lhppp filename\n\nDescription: lhppp processes the data provided by lhp.\n\nfilename: name of the output file generated by lhp that is to be processed ");
        return 1;
    } else {
        fname=argv[1];
    }
    process(fname);
    return 0;
}