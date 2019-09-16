/***************************************************************************
 * Sterile neutrino global fits - utility functions                        *
 ***************************************************************************
 * Author: Joachim Kopp, jkopp@cern.ch                                     *
 ***************************************************************************/
#include <stdio.h>
#include <string.h>
#include "glb_error.h"
#include "nu.h"


// -------------------------------------------------------------------------
int LoadNd(const char *filename, double **buffer, const int n_columns, int *n_rows)
// -------------------------------------------------------------------------
// Read N-dimensional data set from file and store it in buffer. The number
// N == n_columns must be given, while the number of rows, n_rows, is
// returned.  It is the users responsibility to ensure that buffer is large
// enough to hold the whole file.
// -------------------------------------------------------------------------
{
  FILE *f = fopen(filename, "r");
  if (!f)
    return GLBERR_FILE_NOT_FOUND;

  // Read from file
  char line[200];
  *n_rows = 0;
  while (fgets(line, 100, f))
  {
    if (line[0] != '#'  &&  line[strspn(line, " \t")] != '\n')
    {                         // Ignore comments and blank lines
      char *p = line;
      for (int j=0; j < n_columns; j++)
        buffer[*n_rows][j] = strtod(p, &p);
      (*n_rows)++;
    }
  };
  fclose(f);

  if (debug_level > 2)
  {
    printf("Data from file %s:\n", filename);
    for (int i=0; i < *n_rows; i++)
    {
      printf("%5d ", i);
      for (int j=0; j < n_columns; j++)
        printf("%10.5g ", buffer[i][j]);
      printf("\n");
    }
  }

  return GLB_SUCCESS;
}


// -------------------------------------------------------------------------
int LoadNdAlloc(const char *filename, double **buffer, const int n_columns, int *n_rows)
// -------------------------------------------------------------------------
// Read N-dimensional data set from file and store it in buffer. The number
// N == n_columns must be given, while the number of rows, n_rows, is
// returned.  The function allocates memory for buffer as needed.
// -------------------------------------------------------------------------
{
  if (!filename || !n_rows || !buffer)
    return GLBERR_INVALID_ARGS;

  int j;
  int nl;  // Counter for "real" lines, including blank lines and comments
  int buf_length = 100;

  // Allocate buffer
  if (!buffer)
    buffer = new double *[n_columns];
  if  (buffer == NULL)
  {
    fprintf(stderr, "LoadNdAlloc: out of memory while reading file %s.\n", filename);
    return GLBERR_MALLOC_FAILED;
  }
  for (j=0; j < n_columns; j++)
  {
    buffer[j] = new double[buf_length];
    if (buffer[j] == NULL)
    {
      fprintf(stderr, "LoadNdAlloc: out of memory while reading file %s.\n", filename);
      return GLBERR_MALLOC_FAILED;
    }
  }

  FILE *f = fopen(filename, "r");
  if (!f)
  {
    fprintf(stderr, "LoadNdAlloc: Cannot open file %s.\n", filename);
    return GLBERR_FILE_NOT_FOUND;
  }

  // Read from file
  const int max_line = 1024;
  char this_line[max_line];
  *n_rows = nl = 0;
  while (fgets(this_line, max_line, f))
  {
    nl++;
    if (strlen(this_line) > max_line - 2)
    {
      fprintf(stderr, "LoadNdAlloc: Line %d too long in file %s.\n",
              nl, filename);
      return GLBERR_INVALID_FILE_FORMAT;
    }

    if (this_line[0] != '#'  &&  this_line[strspn(this_line, " \t")] != '\n')
    {                               // Ignore comments and blank lines
      if (*n_rows >= buf_length)   // If necessary, resize x-sec buffer
      {
        buf_length *= 2;
        for (j=0; j < n_columns; j++)
        {
          buffer[j] = (double *) realloc(buffer[j], buf_length*sizeof(buffer[j][0]));
          if (buffer[j] == NULL)
          {
            fprintf(stderr, "LoadNdAlloc: out of memory while reading file %s.\n", filename);
            return GLBERR_MALLOC_FAILED;
          }
        }
      }

      char *p = strtok(this_line, " \t,");
      for (j=0; j < n_columns; j++)
      {
        if (!p)
        {
          fprintf(stderr, "LoadNdAlloc: Line %d too short in file %s.\n",
                  nl, filename);
          return GLBERR_INVALID_FILE_FORMAT;
        }
        if (sscanf(p, "%lf", &buffer[j][*n_rows]) != 1)
        {
          fprintf(stderr, "LoadNdAlloc: Line %d invalid in file %s.\n",
                  nl, filename);
          return GLBERR_INVALID_FILE_FORMAT;
        }
        p = strtok(NULL, " \t,");
      }
      (*n_rows)++;
    }
  };

  fclose(f);
  return GLB_SUCCESS;
}

