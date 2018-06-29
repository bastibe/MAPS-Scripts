/* 
   19 Dec 94 - 15 Sept 95 - 12 June 98
   Alain de Cheveigne, LLF, CNRS/Universite Paris 7.
   
   Debugging wrapper for malloc() functions.
   
   */

#define MWRAP_C
#define MWRAP
#include "mwrap.h"

void mexErrMsgTxt(
    const char	*error_msg	
    );


#define GROWSIZE 8192 * 10	/* number of nodes preallocated in
				   node pool.  Each node is 20 bytes. */
				   

/*
   #define NL	fprintf(stderr,"\n");
   #define DOT	fprintf(stderr,".");
   #define DASH	fprintf(stderr,"-");
   #define OK	\
   {fprintf(stderr, "%s<%d> ", __FILE__, __LINE__);\
      fprintf(stderr,"OK\n");}
      */

#define DIE(s)  do {char c[100]; sprintf(c, "%d %s ", mwrap_line, s); mexErrMsgTxt(c);} while(0)
#define DIE1(s,v1)  do {char c[100]; sprintf(c, "%d %s", mwrap_line, s); mexErrMsgTxt(c);} while(0)
#define DIE2(s,v1,v2) do {char c[100]; sprintf(c, "%d %s", mwrap_line, s); mexErrMsgTxt(c);} while(0)
#define DIE3(s,v1,v2,v3) do {char c[100]; sprintf(c, "%d %s", mwrap_line, s); mexErrMsgTxt(c);} while(0)
#define PRINT3(string, v1, v2, v3)  mexErrMsgTxt((string));
	  
	  
/*
#define DIE(s) \
     do {fprintf(stderr, "%s/%d: ", mwrap_file, mwrap_line);\
	   fprintf(stderr, (s)); fprintf(stderr, "\n");  exit(1);} \
     while (0)
#define DIE1(s, value) \
     do {fprintf(stderr, "%s/%d: ", mwrap_file, mwrap_line);\
           fprintf(stderr, (s), (value)); fprintf(stderr, "\n"); exit(1);} \
     while (0)
#define DIE2(s, v1,v2) \
     do {fprintf(stderr, "%s/%d: ", mwrap_file, mwrap_line);\
           fprintf(stderr, s, (v1),(v2)); fprintf(stderr, "\n"); exit(1);} \
     while (0)
#define DIE3(s,v1,v2,v3) \
     do {fprintf(stderr, "%s/%d: ", mwrap_file, mwrap_line);\
	   fprintf(stderr, s, (v1), (v2), (v3)); \
	     fprintf(stderr, "\n"); exit(1);} while (0)
#define PRINT3(string, v1, v2, v3) do {fprintf(stderr,string, v1,v2,v3);} \
     while(0)

		 */

#define NOK(n) mwrap_checknode(n);
     
/* MWRAP_SPLAY() defined as a macro to allow slightly faster splaying, 
   especially when the node is already root or next to root */

#define MWRAP_SPLAY(n) do { \
			mwrap_node father; \
			father = (n)->up; \
			while (father) { \
			  if (!father->up) MWRAP_SINGLE_ROTATE(n); \
			  else MWRAP_DOUBLE_ROTATE(n); \
			  father = (n)->up; \
		        } \
		      } while (0);

#define MWRAP_SINGLE_ROTATE(n) do { \
			mwrap_node p, B; \
		        p = (n)->up; \
			if (p->left == (n)) { \
			  B = (n)->right; (n)->right = p; \
			  p->left = B; if (B) B->up = p; \
			} else { \
			  B = (n)->left; (n)->left = p; \
			  p->right = B; if (B) B->up = p; \
			} \
			p->up = (n); \
			n->up = NULL; \
		      } while (0)

#define MWRAP_DOUBLE_ROTATE(n) mwrap_double_rotate(n)

/* update lowest and highest known stack address */
#define STACK(p) do {\
	       if (mwrap_stack_lo>(char *)(p)) mwrap_stack_lo=(char *)(p); \
	       if (mwrap_stack_hi<(char *)(p)) mwrap_stack_hi=(char *)(p); \
	   } while (0)
     
static void mwrap_register(char *lo, char *hi);
static void mwrap_delete(mwrap_node n);
static mwrap_node mwrap_leftmost(mwrap_node tree);
static void mwrap_checksub(mwrap_node tree, char *low, char *high);
static mwrap_node mwrap_newnode();
static void mwrap_freenode(mwrap_node n) ;
static void mwrap_single_rotate(mwrap_node n);
static void mwrap_double_rotate(mwrap_node n);
static void  mwrap_splay(mwrap_node n);
static void mwrap_check_stack(char *p); /* check if pointer to stack */

static char *heap_lowest; /* lowest address so far */
static char *heap_highest; /* highest address so far +1 */
static mwrap_node tree_lowest; /* lowest address used by tree */
static mwrap_node tree_highest; /* highest address used by tree +1 */
void mwrap_checktree(char *file, long line);

/* global variables */
mwrap_node mwrap_tree = NULL;	/* the tree */
char *mwrap_file;		/* caller source file */
long mwrap_line;		/* line in caller source file */
char *mwrap_stack_lo = (char *) -1; /* lowest known stack address */
char *mwrap_stack_hi = (char *) -1; /* highest known stack address */

/* Mwrap wraps the four basic malloc routines: 
   malloc, calloc , realloc, free.	If you wish it to handle more 
   exotic routines, turn on the following switches. */
     
#define MWRAP_CFREE		/* cfree() */
#define MWRAP_VMALLOC		/* vmalloc() */
#define MWRAP_VFREE		/* vfree() */
/* #define MWRAP_MEMALIGN	/* memalign */
     
/* 
   wrappers for routines of the malloc family 
   */

/* malloc */
DATATYPE *mwrap_malloc(SIZETYPE size, char *file, long line) {
  DATATYPE *p;
  int dummy; STACK(&dummy);

  mwrap_file = file;
  mwrap_line = line;
  p = malloc(size);
  if (!p) DIE1("malloc failed (size=%d)", size);
  /* check for negative size or overflow */
  if ((char *) p > (char *) p+size) 
    DIE1("malloc called with size = %d", size);
  mwrap_register((char *) p, ((char *) p + size));
  return p;
}

/* free */
VOIDTYPE mwrap_free(DATATYPE *p, char *file, long line) {
  int dummy; STACK(&dummy);
  /*fprintf(stderr, "%#x\n", p);  */
  mwrap_file = file;
  mwrap_line = line;
  mwrap_bok(p, file, line);	/* side effect: node is root of mwrap_tree */
  mwrap_delete(mwrap_tree);	/* delete that node */
  free(p);
}

/* realloc */
DATATYPE *mwrap_realloc(DATATYPE *p, SIZETYPE size, char *file, long line) {
  int dummy; STACK(&dummy);
  mwrap_file = file;
  mwrap_line = line;
  mwrap_bok(p, file, line);	/* side effect: node is root of mwrap_tree */
  mwrap_delete(mwrap_tree);	/* delete that node */
  p = realloc(p, size);
  if (!p) DIE1("ralloc failed (size=%d)", size);
  if ((char *) p > (char *) p+size) 
    DIE1("realloc called with size = %d", size);
  mwrap_register((char *) p, (char *) p + size);
  return p;
}

/* calloc */
DATATYPE *mwrap_calloc(SIZETYPE n, SIZETYPE size, char *file, long line) {
  DATATYPE *p;
  int dummy; STACK(&dummy);

  mwrap_file = file;
  mwrap_line = line;
  if (size <= 0) DIE1("calloc called with size = %d", size);
  p = calloc(n, size);
  if (!p) DIE2("calloc failed (n=%d, size=%d)", n, size);
  mwrap_register((char *) p, (char *) p + n*size);
  return p;
}

/* cfree */
#ifdef MWRAP_CFREE
VOIDTYPE mwrap_cfree(DATATYPE *p, char *file, long line) {
  int dummy; STACK(&dummy);
  mwrap_file = file;
  mwrap_line = line;
  mwrap_bok(p, file, line);	/* side effect: node is root of mwrap_tree */
  mwrap_delete(mwrap_tree);
  /* cfree(p); */ free(p);
}
#endif

/* valloc */
#ifdef MWRAP_VALLOC
DATATYPE *mwrap_valloc(SIZETYPE size, char *file, long line) {
  DATATYPE *p;
  int dummy; STACK(&dummy);

  mwrap_file = file;
  mwrap_line = line;
  p = (DATATYPE *) valloc(size);
  if (!p) DIE1("valloc failed (size=%d)", size);
  if ((char *) p > (char *p) p+size) 
    DIE1("valloc called with size = %d", size);
  mwrap_register((char *) p, (char *) p + size);
  return p;
}
#endif

/* vfree */
#ifdef MWRAP_VFREE
VOIDTYPE mwrap_vfree(DATATYPE *p, char *file, long line) {
  int dummy; STACK(&dummy);
  mwrap_file = file;
  mwrap_line = line;
  mwrap_bok(p, file, line);	/* side effect: node is root of mwrap_tree */
  mwrap_delete(mwrap_tree);
  /* vfree(p); */ free(p);
}
#endif

/* memalign */
#ifdef MWRAP_MEMALIGN
DATATYPE *mwrap_memalign(SIZETYPE alignment, SIZETYPE size, char *file, long line) {
  DATATYPE *p;
  int dummy; STACK(&dummy);

  mwrap_file = file;
  mwrap_line = line;

  p = (DATATYPE *) memalign(alignment, size); 
  /* dunno why I need the cast... */
  if (!p) DIE2("memalign failed (alignment=%d, size=%d)", alignment, size);
  if ((char *) p > (char *p) p+size) 
    DIE1("memalign called with size = %d", size);
  mwrap_register((char *) p, (char *) p + size);
  return p;
}
#endif

/* unwrapped versions */

DATATYPE *nowrap_malloc(SIZETYPE size) {
  return malloc(size);
}

/* free */
VOIDTYPE nowrap_free(DATATYPE *p) {
  free(p);
}

/* realloc */
DATATYPE *nowrap_realloc(DATATYPE *p, SIZETYPE size) {
  return realloc(p, size);
}

/* calloc */
DATATYPE *nowrap_calloc(SIZETYPE n, SIZETYPE size) {
  return calloc(n, size);
}

/* cfree */
#ifdef MWRAP_CFREE
VOIDTYPE nowrap_cfree(DATATYPE *p) {
  /* cfree(p); */ free(p);
}
#endif

/* valloc */
#ifdef MWRAP_VALLOC
DATATYPE *nowrap_valloc(SIZETYPE size) {
  return (DATATYPE *) valloc(size);
}
#endif

/* vfree */
#ifdef MWRAP_VFREE
VOIDTYPE nowrap_vfree(DATATYPE *p) {
  /* vfree(p); */ free(p);
}
#endif

/* memalign */
#ifdef MWRAP_MEMALIGN
DATATYPE *nowrap_memalign(SIZETYPE alignment, SIZETYPE size) {
  return (DATATYPE *) memalign(alignment, size); /* why do i need the cast? */
}
#endif

/*     
   routines for handling the binary splay tree (mwrap_tree) that 
   keeps track of limits of allocated blocks
   */

/* register new block in tree */
void mwrap_register(char *lo, char *hi)
{
  mwrap_node n, new;

  /* fprintf(stderr, "%#x %#x %d\n", lo, hi, hi-lo); */
  
  if (lo > hi) 
    DIE2("mwrap_register: low = %#x > hi = %#x (size overflowed?)", lo, hi);
  
  /* alloc and init new node */
  new = mwrap_newnode();
  new->lo = lo;
  new->hi = hi;
  new->right = new->left = NULL;

  /* insert in tree */
  if (mwrap_tree == NULL) {
    new -> up = NULL;
    mwrap_tree = new;
  } else {
    n = mwrap_tree;
    while (n) {
      /* find which side to insert new node */
      if (lo >= n->lo) {	/* must be right */
	if (lo < n->hi) {
	  PRINT3("malloc returned a block (%#x, %#x, size: %d) that overlaps\n", 
		 lo, hi, hi-lo);
	  PRINT3("with previous recorded block (%#x, %#x, size: %d)\n", 
		 n->lo, n->hi, n->hi-n->lo);
	  MWRAPOK();
	  fprintf(stderr, 
		  "malloc heap or mwrap tree may be corrupt!\n"); exit(1);
	}
	if (n->right) {		/* try right subtree */
	  n = n->right;
	} else {		/* put new node here */
	  n->right = new;
	  new->up = n;
	  break;
	}
      } else {			/* must be left */
	if (hi > n->lo) {
	  PRINT3("malloc returned a block (%#x, %#x, size: %d) that overlaps\n", 
		 lo, hi, hi-lo);
	  PRINT3("with previous recorded block (%#x, %#x, size: %d)\n", 
		 n->lo, n->hi, n->hi-n->lo);
	  MWRAPOK();
	  fprintf(stderr, 
		  "malloc heap or mwrap tree may be corrupt!\n"); exit(1);
	}
	if (n->left)		/* try left subtree */
	  n = n->left;
	else {			/* put new node here */
	  n->left = new;
	  new->up = n;
	  break;
	}
      }
    }
    MWRAP_SPLAY(new);		/* bring it to root */
    mwrap_tree = new;
  }

  /* update lowest and highest addresses */
  if (!heap_lowest		/* first call */
      || lo < heap_lowest) heap_lowest = lo;
  if (!heap_highest		/* first call */
      || hi > heap_highest) heap_highest = hi;	
}

/* return leftmost node */
static mwrap_node mwrap_leftmost(mwrap_node tree)
{
  mwrap_node n = tree;
  
  while (n) {
    if (n->left == NULL)
      return n;
    else
      n = n->left;
  }
  DIE("mwrap_leftmost(): tree is NULL");
}

/* delete node from mwrap_tree */
static void mwrap_delete(mwrap_node n) {
  mwrap_node left, right;

  /* bring node to root */
  if (!n) DIE("!");
  MWRAP_SPLAY(n);	
  mwrap_tree = n;

  /* If n has a right subtree:
     - detach it, bring its leftmost node to its root, and make this
     the new root of mwrap_tree,
     - if n has a left subtree, attach it to the (empty) left side of
     mwrap_tree, then free n.
     If n has no right subtree, set mwrap_tree to its left subtree,
     then free n. */

  right = n->right;
  left = n->left;

  if (right) {
    right->up = NULL;		/* to limit splay */
    right = mwrap_leftmost(right);	
    MWRAP_SPLAY(right);			
    mwrap_tree = right;		/* this will be the new root */
    if (mwrap_tree->left) 
      DIE("!!");		/* should be empty */
    if (left) {
      mwrap_tree->left = left; 
      left->up = mwrap_tree;
    }
  } else {
    mwrap_tree = left;
    if (mwrap_tree) mwrap_tree->up = NULL;
  }
  mwrap_freenode(n);
}

/* #define MWRAP_MALLOCNODES */

#ifdef MWRAP_MALLOCNODES	/* call malloc/free for each node */
static void mwrap_freenode(mwrap_node n) 
{
  free (n);
}

static mwrap_node mwrap_newnode()
{
  mwrap_node n;

  n = (mwrap_node) malloc(sizeof(struct mwrap_node_struct));
  if (! tree_lowest		/* first call */
      || n < tree_lowest) tree_lowest = n;
  if (! tree_highest		/* first call */
      || n + sizeof(struct mwrap_node_struct) > tree_highest)
    tree_highest = n + sizeof (struct mwrap_node_struct);
}

#else  /* use our own free node list */

static mwrap_node freelist;

static mwrap_node mwrap_newnode()
{
  mwrap_node n;
  mwrap_node b;
  long i;

  if (!freelist) {
    b = (mwrap_node) malloc(GROWSIZE * sizeof(*b));
    if (!b) DIE("mwrap_newnode(): malloc failed");

    if (! tree_lowest		/* first call */
	|| b < tree_lowest) tree_lowest = b;
    if (! tree_highest		/* first call */
	|| b + GROWSIZE * sizeof(*b) > tree_highest)
      tree_highest = b + GROWSIZE * sizeof (*b);

    freelist = &b[0];
    for (i = 0; i < GROWSIZE-1; i++) {
      b[i].up = &b[i+1];       
    }
    b[GROWSIZE-1].up = NULL;
  }
  n = freelist;
  freelist = n->up;
  return n;
}

static void mwrap_freenode(mwrap_node n)
{
  n->up = freelist;
  freelist = n;
}

#endif


/* splay the tree this node belongs to, bringing the node up to the root */  
/* NOTE: this is replaced by macro MWRAP_SPLAY() */
static void  mwrap_splay(mwrap_node n)
{
  mwrap_node father;

  if (n == NULL) DIE("!");
  father = n->up;
  while(father) {
    if (father->up == NULL) MWRAP_SINGLE_ROTATE(n);
    else MWRAP_DOUBLE_ROTATE(n);
    father = n->up;
  }
}

/* rotate one-segment twig ending on n (certainly at top of tree) */
/* NOTE: this is replaced by macro MWRAP_SINGLE_ROTATE() */
static void  mwrap_single_rotate(mwrap_node n)
{
  mwrap_node p, B;
  
  if (n == NULL) DIE("!");
  p = n->up;

  /* process according to shape of twig */
  if (n->up->left == n) {	/* zig_left */
    B = n->right;
    n->right = p;
    p->left = B; if (B) B->up = p;
  } else {			/* zig right */
    B = n->left;
    n->left = p;
    p->right = B; if (B) B->up = p;
  }

  p->up = n;
  n->up = NULL;			/* because n is now root of the tree */
}

/* rotate twig ending with n */
static void  mwrap_double_rotate(mwrap_node n)
{
  mwrap_node p, g, ggp, B, C;

  if (n == NULL) DIE("!");
  p = n->up;
  g = p->up;
  ggp = g->up;

  /* connect n to rest of tree */
  n->up = ggp;
  if (ggp) {
    if (ggp->left == g) ggp->left = n;
    else ggp->right = n;
  }

  /* process according to shape of twig */
  if (p == g->right) {
    if (n == p->right) {	/* zig_zig_right */
      B = p->left;
      C = n->left;
      n->left = p; p->up = n;
      p->left = g; g->up = p;
      g->right = B; if (B) B->up = g;
      p->right = C; if (C) C->up = p;
    } else {			/* zig_zag_right */
      B = n->left; 
      C = n->right;
      n->left = g; g->up = n;
      n->right = p; p->up = n;
      g->right = B; if (B) B->up = g;
      p->left = C; if (C) C->up = p;
    } 
  } else {			
    if (n == p->left) {		/* zig_zig_left */
      B = n->right;
      C = p->right;
      n->right = p; p->up = n;
      p->right = g; g->up = p;
      p->left = B; if (B) B->up = p;
      g->left = C; if (C) C->up = g;
    } else {			/* zig_zag_left */
      B = n->left; 
      C = n->right;
      n->left = p; p->up = n;
      n->right = g; g->up = n;
      p->right = B; if (B) B->up = p;
      g->left = C; if (C) C->up = g;
    }
  }
}

/* Checking and query routines */

/* check that node address is within plausible range */
static void mwrap_checknode(mwrap_node n)
{
  if (n < tree_lowest) 
    DIE1("checktree: node address %#x is too low to be true", n);
  if (n > tree_highest) 
    DIE1("checktree: node address %#x is too high to be true", n);
}


/* check mwrap_tree for sanity */
void mwrap_checktree(char *file, long line)
{
  mwrap_file = file;
  mwrap_line = line;
  if (!mwrap_tree) return;	/* empty */
  NOK(mwrap_tree);
  if (mwrap_tree) {
    if (mwrap_tree->up) DIE("checktree: mwrap-tree->up is not null");
    mwrap_checksub(mwrap_tree, heap_lowest, heap_highest);
  }
}

/* Recursively check that a subtree has consistent links and backlinks,
   and that block ranges are ordered and do not overlap.  
   At each call the routine is passed two limits, 'low' and 'high'.
   The node's block must be within these limits, and so must blocks in 
   the left and right subtrees.  In addition, blocks in the left
   subtree must be below this node's block, and those in the right
   subtree must be above.  Passing 0 for 'low' or 'high' disables the
   corresponding test (not really necessary...).
   To avoid excessive recursion depth, recursion occurs only when 
   the tree branches. */

static void mwrap_checksub(mwrap_node tree, char *low, char *high)
{
  mwrap_node n;
  char *lo, *hi;

  n = tree;
  while(n) {

    /* check this node */
    NOK(n);
    if (n->lo > n->hi) DIE2("checktree: block limits reversed: %#x > %#x",
			    n->lo, n->hi);
    if (high && n->hi > high) 
      DIE("checktree: block ranges overlap!");
    if (low && n->lo < low) 
      DIE("checktree: block ranges overlap!!");

    /* move along tree until it branches */
    if (!n->left) {
      lo = n->hi;
      if (n->right) {
	NOK(n->right);
	if (n->right->up != n) DIE("checktree: bad link!");
      }
      n = n->right;
    } else if (!n->right) {
      hi = n->lo;
      if (n->left) {
	NOK(n->left);
	if (n->left->up != n) DIE("checktree: bad link!!");
      }
      n = n->left;
    } else {		       
      /* tree branches: recurse on subtrees */
      NOK(n->left);
      if (n->left->up != n) DIE("checktree: bad link!!!");
      mwrap_checksub(n->left, low, n->lo);
      NOK(n->right);
      if (n->right->up != n) DIE("checktree: bad link!!!!");
      mwrap_checksub(n->right, n->hi, high);
      return;
    }
  }
}

/* return base of block containing pointer, NULL if none */
char *mwrap_pbase(char *p)
{
  mwrap_node n = mwrap_tree; 

  while (1) { 
    if (!n) return (NULL);
    if ((char *)(p) < n->lo) 
      n = n->left; 
    else if ((char *)(p) >= n->hi) 
      n = n->right; 
    else break;			/* found */
  } 
  MWRAP_SPLAY(n); 
  mwrap_tree = n; 
  return n->lo;
}

/* return one beyond end of block containing pointer, NULL if none */
char *mwrap_pbump(char *p)
{
  mwrap_node n = mwrap_tree; 

  while (1) { 
    if (!n) return (NULL);
    if ((char *)(p) < n->lo) 
      n = n->left; 
    else if ((char *)(p) >= n->hi) 
      n = n->right; 
    else break;			/* found */
  } 
  MWRAP_SPLAY(n); 
  mwrap_tree = n; 
  return n->hi;
}

/* return size of block containing pointer, -1 if none */
long mwrap_psize(char *p)
{
  mwrap_node n = mwrap_tree;

  while (1) { 
    if (!n) return (-1);
    if ((char *) (p) < n->lo) 
      n = n->left; 
    else if ((char *) (p) >= n->hi) 
      n = n->right; 
    else break;			/* found */
  } 
  MWRAP_SPLAY(n); 
  mwrap_tree = n; 
  return (n->hi - n->lo);
}

/* complain that MWRAP_ON should be defined to use certain macros */
long mwrap_squeal(char *file, long line)
{
  fprintf(stderr, 
	  "%s:%d you must compile with MWRAP_ON defined to use this\n", 
	  file, line);
  exit(1);
}


/* check that a pointer is within an allocated block */
int mwrap_pok(char *p, char *file, long line)
{
  mwrap_node n; 
  
  mwrap_file = file;
  mwrap_line = line;  
  n = mwrap_tree;
  
  /* if (n && p >= n->lo && p < n->hi) return; */

  while (1) { 
    if (!n) { 
      mwrap_checktree(mwrap_file, mwrap_line);
      mwrap_check_stack(p);
      DIE1("%#x is not in allocated block", p);
    }
    if (p < n->lo) n = n->left; 
    else if (p >= n->hi) n = n->right; 
    else break; 
  }
  MWRAP_SPLAY(n);		/* side effect: queried block is now root */
  mwrap_tree = n; 
  return 1;
}

/* check that pointer is base of allocated block */
int mwrap_bok(char *b, char *file, long line)
{ 
  mwrap_file = file;
  mwrap_line = line;
  mwrap_pok(b, file, line);	/* side effect: block is root */

  if (b != mwrap_tree->lo) { 
    mwrap_checktree(mwrap_file, mwrap_line); /* tree damaged ? */
    mwrap_check_stack(b);	/* might point to stack? */
    DIE1("%#x is not the base of an allocated block", b ); 
  }
  return 1;
}

/* check that pointers are both in the same block */
int mwrap_pqok(char *p, char *q, char *file, long line)
{ 
  mwrap_file = file;
  mwrap_line = line;
  mwrap_pok(p, file, line);	/* side effect: block is root */
  
  if (q < mwrap_tree->lo || q >= mwrap_tree->hi) { 
    mwrap_checktree(file, line); /* tree damaged ? */
    mwrap_pok(q, file, line);	/* is q at least in a block ? */
    DIE2("%#x and %#x are not in same block",  p, q); 
  } 
} 

/* check that pointer p is in block of base b */
int mwrap_pbok(char *p, char *b, char *file, long line)
{ 
  mwrap_file = file;
  mwrap_line = line;
  mwrap_pok(p, file, line);	/* side effect: block is root of tree */
  
  if (b != mwrap_tree->lo) { 
    char *a = mwrap_tree->lo;
    mwrap_checktree(file, line); /* tree damaged ? */
    mwrap_pok(b, file, line);	/* is it at least in a block ? */
    DIE3("%#x != base (%#x) of block containing %#x", 
	 b, a, p); 
  } 
  return 1;
}

/* routines to register or remove blocks explicitly */
void mwrap_checkin(char *lo, char *hi, char *file, long line) {
  mwrap_node n, new;

  mwrap_file = file;
  mwrap_line = line;

  /* This code does the work of mwrap_register, plus some checks */

  /* We can't trust the lo-hi interval not to overlap with 
     existing blocks, so we check here (actually we can't 
     trust it not to overlap with future blocks, but we
     can't do much about that...) */

  if (lo > hi) 
    DIE2("mwrap_checkin: low = %#x > hi = %#x", lo, hi);
  
  /* alloc and init new node */
  new = mwrap_newnode();
  new->lo = lo;
  new->hi = hi;
  new->right = new->left = NULL;

  /* insert in tree */
  if (mwrap_tree == NULL) {
    new -> up = NULL;
    mwrap_tree = new;
  } else {
    n = mwrap_tree;
    while (n) {
      /* find which side to insert new node */
      if (lo >= n->lo) {	/* must be right */
	if (lo < n->hi) 
	  DIE1("mwrap_checkin(): trying to put lo (%#x) in existing block", lo);
	if (n->right) {		/* try right subtree */
	  n = n->right;
	} else {		/* put new node here */
	  n->right = new;
	  new->up = n;
	  if (lo < n->hi)
	    DIE1("mwrap_checkin(): lo (%#x) overlaps with previous block", lo);
	  break;
	}
      } else {			/* must be left */
	if (hi > n->lo) 
	  DIE1("mwrap_checkin(): trying to put hi (%#x) in existing block", hi);
	if (n->left)		/* try left subtree */
	  n = n->left;
	else {			/* put new node here */
	  n->left = new;
	  new->up = n;
	  if (hi > n->lo)
	    DIE1("mwrap_checkin(): hi (%#x) overlaps with next block", hi);
	  break;
	}
      }
    }
    MWRAP_SPLAY(new);		/* bring it to root */
    mwrap_tree = new;
  }

  /* update lowest and highest addresses */
  if (!heap_lowest		/* first call */
      || lo < heap_lowest) heap_lowest = lo;
  if (!heap_highest		/* first call */
      || hi > heap_highest) heap_highest = hi;	
}

void mwrap_checkout(char *lo, char *file, long line) {
  mwrap_file = file;
  mwrap_line = line;
  mwrap_bok(lo, file, line);	/* side effect: node is root of mwrap_tree */
  mwrap_delete(mwrap_tree);	/* delete that node */
}

/* get pointer near top of stack, check if pointer is in stack range */
static void mwrap_check_stack(char *p)
{
  int dummy; STACK(&dummy);

  if (p >= mwrap_stack_lo && p <= mwrap_stack_hi)
    fprintf(stderr, "(could be on stack) ");
}
