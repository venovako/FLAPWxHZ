#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  if (argc != 2) {
    (void)fprintf(stderr, "%s binary_file_name\n", argv[0]);
    return EXIT_FAILURE;
  }

  FILE *const f = fopen(argv[1], "r");
  if (!f) {
    (void)fprintf(stderr, "%s cannot be opened for reading\n", argv[1]);
    return EXIT_FAILURE;
  }

  for (double d; fread(&d, sizeof(d), 1u, f); fprintf(stdout, "%#26.17e\n", d));

  (void)fclose(f);
  return EXIT_SUCCESS;
}
