#define KTUP 12
#define MASK 16777215   // should be (4^KTUP)-1

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

struct feature {
  unsigned feature_index;
  unsigned fragment_index;
  unsigned position;
  unsigned shift;
};

extern int NFEATURES;
extern void load_features (char *);
extern void* get_iterator_data ();
extern void free_iterator_data (void *v);
extern int
iterate (const char *s,
         unsigned size, unsigned idx,
         void *v, struct feature *f);

struct feature_desc {
  unsigned feature_index;
  unsigned fragment_index;
  unsigned mask;
  unsigned seq;
  unsigned shift;
};

struct site_ext { unsigned short on; };

struct feature_desc *F;
int M_zero = -1;
unsigned F_index[MASK+1];

int NFEATURES = 0;
#define MAX_NFEATURES 1024*1024

void
load_features (char *blastdata)
{
  char buf[1024];
  FILE *fp;
  unsigned n;

  fp = fopen (blastdata, "r");
  if (fgets (buf, 1024, fp) == NULL) {
    assert ("bad index");
    return;
  }
  NFEATURES = atoi (buf);
  assert (NFEATURES < MAX_NFEATURES);
  F = (struct feature_desc*)
    malloc (sizeof (struct feature_desc)*(NFEATURES));
  bzero(F_index,sizeof(F_index));

  n = 0;
  while (fgets (buf, 1024, fp) && n < NFEATURES) {
    char word[1024];
    unsigned i, j, nw;
    nw = 0;
    for (i=0,j=0; i<strlen (buf); i++) {
      if (buf[i] != ',')
        word[j++] = buf[i]; 
      else {
        int v;
        word[j] = '\0';
        v = atoi (word);
	    if (nw == 0)
	      F[n].feature_index = v;
	    else if (nw == 1)
	      F[n].fragment_index = v;
	    else if (nw == 2)
	      F[n].mask = v;
	    else if (nw == 3)
	      F[n].seq = v;
	    else if (nw == 4)
	      F[n].shift = v;
        nw++;
	    j=0;
      }
    }
	if (F[n].mask != 0 && M_zero < 0)
	  M_zero = n;
	else if (F_index[F[n].seq] == 0)
	  F_index[F[n].seq]=n+1;
	n++;
  }

  fclose (fp);
}

void* get_iterator_data ()
{
  unsigned n;
  void *x;
  unsigned *y;
  unsigned i;
  struct site_ext *v;
  n = sizeof (struct site_ext) * (NFEATURES);
  n += sizeof (unsigned);
  x = (void *)malloc (n);
  y = (unsigned*)x;
  *y = 0;
  v = (struct site_ext *)(x+sizeof(unsigned));
  for (i=0; i<NFEATURES; i++) { v [i].on = 0; }
  return x;
}

void free_iterator_data (void *v) { free (v); }
  
static int on_list[MAX_NFEATURES];
unsigned int on_list_cur = 0;

int
iterate (const char *s, unsigned size, unsigned idx,
         void *v, struct feature *f)
{
  unsigned i,oi;
  unsigned *p = (unsigned *)v;
  register unsigned sval = (unsigned) (*p);
  register unsigned bidx = 0;
  struct site_ext *x = (struct site_ext *) (v+sizeof(unsigned));
  if (idx >= KTUP-1) bidx = idx-(KTUP-1);
  while (idx < size) {
    for (oi=0; oi<on_list_cur; oi++) {
      if (on_list[oi] == -1) continue;
      i = on_list[oi];
	  on_list[oi] = -1;
      if (x [i].on) {
        f->feature_index = F[i].feature_index;
        f->fragment_index = F[i].fragment_index;
        f->position = bidx-1;
        f->shift = F[i].shift;
        x [i].on = 0;
        *p = sval;
        return idx;
      }
    }
    on_list_cur = 0;

    if (s [idx] == 'A' || s [idx] == 'a')
      sval = (sval << 2) + 0;
    else if (s [idx] == 'G' || s [idx] == 'g')
      sval = (sval << 2) + 1;
    else if (s [idx] == 'C' || s [idx] == 'c')
      sval = (sval << 2) + 2;
    else if (s [idx] == 'T' || s [idx] == 't')
      sval = (sval << 2) + 3;
    if (idx < KTUP-1) {
      idx++;
      continue;
    }
    sval = sval & MASK;

#if 1
    if (F_index[sval]>0) {
      for (i=F_index[sval]-1;
	       i<NFEATURES && F[i].seq == sval && F[i].mask == 0; i++) {
        if (x [i].on == 0) {
          x [i].on = 1;
		  on_list[on_list_cur] = i;
		  on_list_cur++;
        }
      }
	}
	if (M_zero >= 0) {
	  for (i=M_zero; i<NFEATURES; i++) {
        if (x [i].on == 0) {
          if ((sval & F[i].mask) == F[i].seq) {
            x [i].on = 1;
		    on_list[on_list_cur] = i;
		    on_list_cur++;
          }
        }
	  }
	}
#else
    for (i=0; i<NFEATURES; i++) {
      if (x [i].on == 0) {
        if ((F[i].mask == 0 && sval == F[i].seq) ||
	        (F[i].mask && (sval & F[i].mask) == F[i].seq)) {
          x [i].on = 1;
		  on_list[on_list_cur] = i;
		  on_list_cur++;
        }
      }
    }
#endif
    idx++;
    bidx++;
  }
  for (oi=0; oi<on_list_cur; oi++) {
    if (on_list[oi] == -1) continue;
    i = on_list[oi];
	on_list[oi] = -1;
    if (x [i].on) {
      f->feature_index = F[i].feature_index;
      f->fragment_index = F[i].fragment_index;
      f->position = bidx-1;
      f->shift = F[i].shift;
      x [i].on = 0;
      *p = sval;
      return idx;
    }
  }
  on_list_cur = 0;

  return -1;
}

#define MAXSEQLEN 600000

void
get_frags (char *file, FILE *fp)
{
  int c;
  int n = 0;
  int n2 = 0;
  int idx = 0;
  char s[MAXSEQLEN+1];
  char s2[2*MAXSEQLEN+1];
  void *v;

  while ((c = fgetc (fp)) != EOF && n < MAXSEQLEN) {
    // XXX: we don't handle wildcards yet
    if (toupper (c) == '*' || toupper (c) == 'N')
      c = 'A';
    else if (c == ' ' || c == '\n')
      continue;
    else if (toupper (c) != 'A' && toupper (c) != 'G' &&
	         toupper (c) != 'C' && toupper (c) != 'T') {
      c = 'A';
      // fprintf (stderr, "%s cut off at pos %d (%d)\n", file, n, c);
      // break;
    }
    s[n] = toupper (c); 
    n++;
  }
  s[n] = '\0';
  if (c != EOF && n >= MAXSEQLEN)
    fprintf (stderr, "%s cut off at %d (max)\n", file, n);

  printf ("======%d %s\n", n, file);

  s2[0] = '\0';
  strcat (s2, s);
  strcat (s2, s);
  n2 = n+n;

  v = get_iterator_data ();
  while (idx >= 0) {
    struct feature f;
    idx = iterate (s2, n2, idx, v, &f);
    if (idx >= 0 && f.position<n2)
      printf ("%d %d %d %d\n",
	      f.feature_index, f.fragment_index, f.position+1, f.shift);
  }
  free_iterator_data (v);
}

int
main (int argc, char *argv[])
{
  unsigned i;
  if (argc < 3) {
    fprintf (stderr, "usage: %s <blastdata> <sequence files...>\n", argv [0]);
    return 1;
  }
  load_features (argv[1]);
  for (i=2; i<argc; i++) {
    FILE *fp = fopen (argv[i], "r");
    if (fp == NULL)
      continue;
    get_frags (argv[i], fp);
    fclose (fp);
  }
  return 0;
}

