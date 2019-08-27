#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

enum { MAXL = 52, MAXC = 256};

FILE *xfopen (const char *fn, const char *mode);
int badmode (const char *s);
int xfclose (FILE *fp);
int xfexists (char *fn);
char *fnwoext (char *nm, char *fn);

int main (int argc, char **argv) {

    int pid, status;
    float f[MAXL] = {0.0};
    float x_line[MAXL] = {0.0};
    char *fn = argc > 1 ? argv[1] : "gnuplot.dat";
    char fnbase[MAXC] = "", fnplt[MAXC] = "";
    size_t i;
    FILE *fp = NULL;
    double x =  -2.5;

    for (i = 0; i < MAXL; i++){      /* fill array of values   */
        f[i] = x * x * x * x - 5 * x * x - x - 3;   /* x^3 - x^2. no overflow */
        x_line[i] = x;
        x+= 0.1;
    }

    fp = xfopen (fn, "w");      /* open output file */

    for (i = 0; i < MAXL; i++)  /* write values to file */
        fprintf (fp, "%10.2f %10.2f\n", x_line[i], f[i]);

    xfclose (fp);   /* close output file */

    /* create 'plot' file 'fn.plt' */
    strcpy (fnplt, fnwoext (fnbase, fn));
    strcat (fnplt, ".plt");
    if (!xfexists (fnplt)) {
        xfopen (fnplt, "w");
        fprintf (fp, "set xlabel 'x'\n"
                    "set ylabel 'f(x) = x^3 - x^2'\n"
                    "set title 'Function Plot of f(x) = x^3 - x^2'\n"
                    "set grid\n"
                    "set style data lines\n"
                    "plot \"%s\" using 1:2 lw 3 linecolor rgb \"blue\"\n"
                    "quit\n", fn);
        xfclose (fp);
    }

    /* fill arguments array for execvp */
    char *args[] = { "gnuplot", "-p", fnplt, NULL };

    if ((pid = (fork())) < 0) { /* fork plot process */
        fprintf (stderr, "fork() error: fork failed.\n");
        return 1;
    }
    else if (pid == 0) {    /* plot from child process */
        if (execvp (*args, args) == -1) {
            fprintf (stderr, "execvp() error: returned error.\n");
            _exit (EXIT_FAILURE);
        }
    }

    waitpid (pid, &status, 0);  /* wait for plot completion (not req'd) */

    return 0;
}

/** fopen with error checking - short version */
FILE *xfopen (const char *fn, const char *mode)
{
    if (!fn || !mode || badmode (mode)) {
        fprintf (stderr, "xfopen() error: invalid parameter.\n");
        exit (EXIT_FAILURE);
    }
    FILE *fp = fopen (fn, mode);

    if (!fp) {
        fprintf (stderr, "xfopen() error: file open failed '%s'.\n", fn);
        exit (EXIT_FAILURE);
    }

    return fp;
}

/** validate file mode 's' is "rwa+b" */
int badmode (const char *s)
{
    const char *modes = "rwa+b";

    for (; *s; s++) {
        const char *m = modes;
        int valid = 0;
        while (*m) if (*s == *m++) { valid = 1; break; }
        if (!valid) return *s;
    }
    return 0;
}

/** file close with error check */
int xfclose (FILE *fp)
{
    if (fclose (fp)) {
        fprintf (stderr, "xfclose() error: nonzero return on fclose.\n");
        return 1;
    }
    return 0;
}

/** check if file 'fn' already exists */
int xfexists (char *fn)
{
    /* if access return is not -1 file exists */
    if (access (fn, F_OK ) != -1 )
        return 1;

    return 0;
}

/** isolate filename, without path or extension */
char *fnwoext (char *nm, char *fn)
{
    char *p  = NULL, *ep = NULL;
    char fnm[MAXC] = "";

    if (!fn) return NULL;
    strcpy (fnm, fn);
    if ((p = strrchr (fnm, '/')))
        p++;
    else
        p = fnm;

    if ((ep = strrchr (p, '.'))) {
        *ep = 0;
        strcpy (nm, p);
        *ep = '.';
    } else
        strcpy (nm, p);

    return nm;
}
