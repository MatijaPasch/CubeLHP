#include <stdlib.h>
#include <stdio.h>

typedef unsigned long ulong;
typedef unsigned char uchar;
typedef struct vertstruct vert;
typedef struct hedgstruct hedg;
typedef struct cubestruct cube;
typedef struct nodestruct node;
typedef struct loosepathstruct loosepath;
typedef struct loosepathsstruct loosepaths;

struct vertstruct {
    uchar *label;
    hedg **e;
    uchar covered;
};

struct hedgstruct {
    vert **v;
    uchar free;
};

struct cubestruct {
    uchar k;
    uchar d;
    vert *v;
    hedg *e;
    ulong vcnt;
    ulong ecnt;
    ulong lhplen;
};

struct nodestruct {
    hedg *e;
    vert *v;
    uchar d;
    uchar *vals;
    ulong blockedcnt;
    hedg** blocked;
};

struct loosepathstruct {
    ulong len;
    hedg **e;
    vert **v;
    loosepath *next;
};

struct loosepathsstruct {
    FILE *f;
    ulong cnt;
    ulong *lcnt;
    loosepath *first;
    loosepath *last;
    ulong lmin;
    ulong lmax;
};

ulong powuc(uchar a, uchar b) {
    ulong res;
    uchar i;
    res=1;
    for(i=0;i<b;i++) {
        res=res*a;
    }
    return res;
}

ulong powuluc(ulong a, uchar b) {
    ulong res;
    uchar i;
    res=1;
    for(i=0;i<b;i++) {
        res=res*a;
    }
    return res;
}

ulong factorial(uchar a) {
    ulong res;
    uchar b;
    res=1;
    for(b=2;b<=a;b++) {
        res*=b;
    }
    return res;
}

cube *initcube(uchar k, uchar d) {
    ulong i,ind,ind2,p;
    uchar coord,coord2, val, *elabel;
    cube *c=malloc(sizeof(cube));
    c->k=k;
    c->d=d;
    c->vcnt=powuc(k,d);
    c->ecnt=d*powuc(k,d-1);
    c->lhplen=(c->vcnt-1)/(k-1);
    c->v=malloc(c->vcnt*sizeof(vert));
    c->e=malloc(c->ecnt*sizeof(hedg));
    for(i=0;i<c->vcnt;i++) {
        c->v[i].covered=0;
        c->v[i].label=malloc(d*sizeof(uchar));
        ind=i;
        for(coord=0;coord<d;coord++) {
            ind2=ind/k;
            c->v[i].label[coord]=ind-k*ind2;
            ind=ind2;
        }
        c->v[i].e=malloc(d*sizeof(hedg *));
        for(coord=0;coord<d;coord++) {
            ind=0;
            p=1;
            for(coord2=0;coord2<coord;coord2++) {
                ind=ind+c->v[i].label[coord2]*p;
                p=p*k;
            }
            for(coord2=coord+1;coord2<d;coord2++) {
                ind=ind+c->v[i].label[coord2]*p;
                p=p*k;
            }
            ind=ind+coord*p;
            c->v[i].e[coord]=(c->e)+ind;
        }

    }
    elabel=malloc(d*sizeof(uchar));
    for(i=0;i<c->ecnt;i++) {
        c->e[i].free=1;
        c->e[i].v=malloc(k*sizeof(vert *));
        ind=i;
        for(coord=0;coord<d-1;coord++) {
            ind2=ind/k;
            elabel[coord]=ind-k*ind2;
            ind=ind2;
        }
        elabel[d-1]=ind;
        for(val=0;val<k;val++) {
            ind=0;
            p=1;
            for(coord=0;coord<elabel[d-1];coord++) {
                ind=ind+elabel[coord]*p;
                p=p*k;
            }
            ind=ind+val*p;
            p=p*k;
            for(coord=elabel[d-1]+1;coord<d;coord++) {
                ind=ind+elabel[coord-1]*p;
                p=p*k;
            }
            c->e[i].v[val]=(c->v)+ind;
        }
    }
    free(elabel);
    return c;
}

void freecube(cube *c) {
    ulong i;
    for(i=0;i<c->ecnt;i++) {
        free(c->e[i].v);
    }
    for(i=0;i<c->vcnt;i++) {
        free(c->v[i].label);
        free(c->v[i].e);
    }
    free(c->e);
    free(c->v);
    free(c);
}

node *initnode(cube *c) {
    uchar i;
    node *n=malloc(sizeof(node));
    n->vals=malloc(c->d*sizeof(uchar));
    n->blockedcnt=0;
    n->blocked=malloc((c->k-1)*(c->d-1)*sizeof(hedg *));
    return n;
}

void freenode(node *n) {
    free(n->blocked);
    free(n->vals);
    free(n);
}

node **initnodes(cube *c,loosepaths *lps) {
    ulong i;
    // we will consider all edgs from 1st to last (even though first two are fixed)
    node **ns=malloc(lps->lmax*sizeof(node*));
    for(i=0;i<lps->lmax;i++) ns[i]=initnode(c);
    ns[0]->e=c->e; //1st edg, containing origin & variation of 1st coord
    ns[0]->v=c->v+1; //1 on 1st, 0 else
    ns[0]->d=1; //only one coord touched so far
    ns[0]->vals[0]=2; //on coord 0 we have taken 0 & 1
    return ns;
}

void freenodes(loosepaths *lps, node **ns) {
    ulong i;
    for(i=0;i<lps->lmax;i++) freenode(ns[i]);
    free(ns);
}

void enternode(cube *c, node *n) {
    uchar val,coord;
    vert *v;
    // block the edg we just took
    n->e->free=0;
    for(val=0;val<c->k;val++) {
        // for all verts in new edg but new endpt cover them (incl old endpt which now has deg 2)
        if(n->e->v[val]!=n->v) {
            v=n->e->v[val];
            v->covered=1;
            // also, block all free inc edgs (notice that current edg is blocked in 1st line, so at most d-1 blocked for each of k-1 non-endpts)
            for(coord=0;coord<c->d;coord++) {
                if(v->e[coord]->free) {
                    v->e[coord]->free=0;
                    n->blocked[n->blockedcnt]=v->e[coord];
                    n->blockedcnt+=1;
                }
            }
        }
    }
}

void leavenode(cube *c, node *n) {
    ulong i;
    for(i=0;i<n->blockedcnt;i++) n->blocked[i]->free=1;
    n->blockedcnt=0;
    n->e->free=1;
    // covered all except endpt, so should uncover all except endpt - but: endpt was uncovered before we took edg & after we took edg, so it's not covered anyway
    for(i=0;i<c->k;i++) n->e->v[i]->covered=0;
}

loosepath *initloosepath(ulong len) {
    loosepath *lp=malloc(sizeof(loosepath));
    lp->len=len;
    if(len>0) {
        lp->v=malloc(len*sizeof(vert *));
        lp->e=malloc(len*sizeof(hedg *));
    } else {
        lp->v=NULL;
        lp->e=NULL;
    }
    lp->next=NULL;
    return lp;
}

void freeloosepath(loosepath *lp) {
    if(lp->len>0) {
        free(lp->v);
        free(lp->e);
    }
    free(lp);
}

void lp_add(loosepaths *lps, ulong len, node **ns) {
    ulong i;
    loosepath *lp=initloosepath(len);
    for(i=0;i<len;i++) {
        lp->e[i]=ns[i]->e;
        lp->v[i]=ns[i]->v;
    }
    lps->last->next=lp;
    lps->last=lp;
    lps->lcnt[len-lps->lmin]+=1;
    lps->cnt+=1;

}


void updatefile(cube *c,loosepaths *lps) {
    ulong i,j;
    loosepath *lp;
    vert *stpt,*v;
    uchar coord;
    lp=lps->last;
    fprintf(lps->f,"%lu,",lp->len);
    stpt=c->v;
    for(coord=0;coord<c->d;coord++) fprintf(lps->f,"%u",stpt->label[coord]);
    for(i=0;i<lp->len;i++) {
        for(j=0;j<c->k;j++) {
            v=lp->e[i]->v[j];
            if(v!=stpt && v!=lp->v[i]) {
                fprintf(lps->f,",");
                for(coord=0;coord<c->d;coord++) fprintf(lps->f,"%u",v->label[coord]);
            }
        }
        fprintf(lps->f,",");
        for(coord=0;coord<c->d;coord++) fprintf(lps->f,"%u",lp->v[i]->label[coord]);
        stpt=lp->v[i];
    }
    fprintf(lps->f,"\n");
    fflush(lps->f);
}

void storelhp(cube *c, node **ns, loosepaths *lps,ulong h) {
        lp_add(lps,h,ns);
        printf("Total number of loose paths: %u\t\r",lps->cnt);
        updatefile(c,lps);
}

void dfs(cube *c, node **ns, loosepaths *lps, ulong h) {
    node *n;
    uchar coord,val;
    if(h<c->lhplen) {
        n=ns[h-1];
        enternode(c,n);
        // these are the correct values for the next nodes if we visit old coords & old vals
        ns[h]->d=n->d;
        for(coord=0;coord<n->d;coord++) ns[h]->vals[coord]=n->vals[coord];
        // we start by visiting old coords
        for(coord=0;coord<n->d;coord++) {
            // if the hedg along this coord is free, we take all verts in hedg except current as endpts (oldvals & new val)
            if(n->v->e[coord]->free) {
                // this will be the edg forall endpts
                ns[h]->e=n->v->e[coord];
                for(val=0;val<n->vals[coord];val++) {
                    // not val of current vert n->v (on this coord)
                    if(val!=n->v->label[coord]) {
                        // set vert for oldval
                        ns[h]->v=ns[h]->e->v[val];
                        // d & vals & e & v all set, ready to go
                        dfs(c,ns,lps,h+1);
                    }
                }
                // time to try a new val, if there's one left
                if(n->vals[coord]<c->k) {
                    ns[h]->vals[coord]+=1; //since we take new val, it will be an oldval from there on
                    ns[h]->v=ns[h]->e->v[n->vals[coord]]; //the vertex for the (current) new val
                    // d & vals & e & v all set, ready to go
                    dfs(c,ns,lps,h+1);
                    ns[h]->vals[coord]-=1; //we're done with the branch where we took the new val, back to default for upcoming branches
                }
            }
        }
        // time to try a new coord, if there's one left
        if(n->d<c->d) {
            ns[h]->d+=1; //since we take a new coord, it will be an oldcoord from there on
            ns[h]->vals[n->d]=2; //change to new coord means curval is zero, have to change to new val, which is one, will be oldval, so newval=2
            ns[h]->e=n->v->e[n->d]; //new edg is the variation of new coord
            ns[h]->v=ns[h]->e->v[1]; //new vert is vert in new edg with val 1 in free coord (which is new coord)
            dfs(c,ns,lps,h+1);
            // no need to clean up, we move to lower level anyway
        }
        leavenode(c,n);
    }
    if(h>=lps->lmin&&h<=lps->lmax) storelhp(c,ns,lps,h);
}

FILE *initfile(cube *c, loosepaths *lps) {
    char fname[256];
    ulong x,i;
    FILE *f;
    sprintf(fname,"lhps_k%u_d%u_dmin%lu_dmax%lu.csv",c->k,c->d,c->lhplen-lps->lmin,c->lhplen-lps->lmax);
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
    fprintf(f,"\nRESULTS\n\nk,%u\nd,%u\n#vertices,%lu\n#edges,%lu\nm,%lu\nm_min,%lu\nm_max,%lu\n",c->k,c->d,c->vcnt,c->ecnt,c->lhplen,lps->lmin,lps->lmax);
    x=factorial(c->d)*powuluc(factorial(c->k-1),c->d);
    fprintf(f,"Type 2 LPs per Type 1 LP,%lu\n",x);
    fprintf(f,"Type 3 LPs per Type 2 LP,%lu\n",c->vcnt);
    x=x*c->vcnt;
    fprintf(f,"Type 3 LPs per Type 1 LP,%lu\n\n",x);
    fprintf(f,"TYPE 1 LPS (displayed as ordered LPs)\n\nLength");
    for(i=0;i<c->vcnt;i++) fprintf(f,",v[%lu]",i);
    fprintf(f,"\n");
    fflush(f);
    return f;
}

loosepaths *initloosepaths(cube *c,ulong lmin,ulong lmax) {
    ulong i;
    loosepaths *lps=malloc(sizeof(loosepaths));
    lps->cnt=0;
    lps->first=initloosepath(0);
    lps->last=lps->first;
    lps->lmin=lmin;
    lps->lmax=lmax;
    lps->lcnt=malloc((lmax-lmin+1)*sizeof(ulong));
    for(i=0;i<=lmax-lmin;i++) lps->lcnt[i]=0;
    lps->f=initfile(c,lps);
    return lps;
}

void finishfile(cube *c,loosepaths *lps) {
    ulong cnt,fac12,fac13,t1,t2,t3;
    fprintf(lps->f,"\nNUMBERS OF DIRECTED LPS (for undirected LPs divide by 2)\n\n");
    fprintf(lps->f,"Length,Type 1,Type 2,Type 3\n");
    fac12=factorial(c->d)*powuluc(factorial(c->k-1),c->d);
    fac13=fac12*c->vcnt;
    for(cnt=0;cnt<=lps->lmax-lps->lmin;cnt++) {
        t1=lps->lcnt[cnt];
        t2=fac12*t1;
        t3=fac13*t1;
        fprintf(lps->f,"%lu,%lu,%lu,%lu\n",lps->lmin+cnt,t1,t2,t3);
    }
    fprintf(lps->f,"Total,%lu,%lu,%lu\n",lps->cnt,lps->cnt*fac12,lps->cnt*fac13);
    fclose(lps->f);
}

void freeloosepaths(cube *c,loosepaths *lps) {
    ulong i;
    loosepath *cur, *next;
    finishfile(c,lps);
    cur=lps->first;
    for(i=0;i<=lps->cnt;i++) {
        next=cur->next;
        freeloosepath(cur);
        cur=next;
    }
    free(lps->lcnt);
    free(lps);
}

void enumeratelhps(cube *c,loosepaths *lps) {
    node **ns=initnodes(c,lps);
    printf("Total number of loose paths: 0\t\r");
    dfs(c,ns,lps,1);
    freenodes(lps,ns);
}

void testcube(uchar k,uchar d) {
    ulong i;
    uchar i2,i3,i4;
    cube *c=initcube(k,d);
    for(i=0;i<c->vcnt;i++) {
        for(i2=0;i2<d;i2++) printf("%u",c->v[i].label[i2]);
        printf("\n");
        for(i2=0;i2<d;i2++) {
            for(i3=0;i3<k;i3++) {
                for(i4=0;i4<d;i4++) printf("%u",c->v[i].e[i2]->v[i3]->label[i4]);
                if(i3<k-1) printf(",");
            }
            printf("\n");
        }
        printf("\n");
    }
    for(i=0;i<c->ecnt;i++) {
        for(i2=0;i2<k;i2++) {
            for(i3=0;i3<d;i3++) printf("%u",c->e[i].v[i2]->label[i3]);
            if(i2<k-1) printf(",");
        }
        printf("\n");
    }
    freecube(c);
}

int main(int argc, char *argv[]) {
    uchar k;
    uchar d;
    ulong dmin,dmax;
    uchar debug=0;
    if(debug) {
        k=2;
        d=4;
        dmin=1;
        dmax=0;
    } else if(argc<3 || argc>5) {
        printf("\nUsage: lhp k d dmin dmax\n\nlhp computes all loose paths (LPs) of a length (=hyperedge count) in the given interval in the cube hypergraph.\n\nk: hyperedge size (=coordinate domain size)\nd: cube hypergraph dimension (=vertex degree)\ndmin: LP length must be at least m-dmin, where m is the loose Hamilton path length (optional, defaults to 0)\ndmax: loose path length must be at most m-dmax (optional, defaults to 0)\n");
        return 1;
    } else {
        sscanf(argv[1],"%hhu",&k);
        sscanf(argv[2],"%hhu",&d);
    }
    cube *c=initcube(k,d);
    dmin=0;dmax=0;
    if(argc>3) sscanf(argv[3],"%lu",&dmin);
    if(dmin>=c->lhplen) {
        printf("WARNING: min has to be smaller than %lu, taking default 0",c->lhplen);
        dmin=0;
    }
    if(argc>4) sscanf(argv[4],"%lu",&dmax);
    if(dmax>dmin) {
        printf("WARNING: max has to be at most %lu, taking default 0",dmin);
        dmax=0;
    }
    loosepaths *lps=initloosepaths(c,c->lhplen-dmin,c->lhplen-dmax);
    enumeratelhps(c,lps);
    freeloosepaths(c,lps);
    freecube(c);
    return 0;
}
