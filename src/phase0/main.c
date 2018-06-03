int main(int argc, char* argv[])
{
  if (argc != 3) {
    (void)fprintf(stderr, "%s input.bin output.bin\n", *argv);
    return EXIT_FAILURE;
  }

  unsigned L = 0u;
  unsigned a = 0u;
  unsigned G = 0u;
  double complex *A = (double complex*)NULL;
  double complex *B = (double complex*)NULL;
  unsigned *T_sizes = (unsigned*)NULL;
  double complex **T = (double complex**)NULL;
  unsigned *lmaxs = (unsigned*)NULL;
  unsigned max_lmax = 0u;
  double *u_norms = (double*)NULL;

  long ret = load_data(argv[1], &L, &a, &G, &A, &B, &T_sizes, &T, &lmaxs, &max_lmax, &u_norms);
#ifndef NDEBUG
  (void)fprintf(stdout, "load_data = %ld\n", ret);
  (void)fflush(stdout);
#endif // !NDEBUG
  if (ret < 0)
    return EXIT_FAILURE;

  ret = split_data(argv[2], L, a, G, max_lmax, T_sizes, A, B, T, lmaxs, u_norms);
#ifndef NDEBUG
  (void)fprintf(stdout, "split_data = %ld\n", ret);
  (void)fflush(stdout);
#endif // !NDEBUG
  return ((ret < 0) ? EXIT_FAILURE : EXIT_SUCCESS);
}
