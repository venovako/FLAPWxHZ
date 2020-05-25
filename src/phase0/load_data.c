static long load_data
(const char *const file,
 unsigned *const L,
 unsigned *const a,
 unsigned *const G,
 double complex **const A,
 double complex **const B,
 unsigned **const T_sizes,
 double complex ***const T,
 unsigned **const lmaxs,
 unsigned *const max_lmax,
 double **const u_norms)
{
  FILE *const fp = fopen(file, "rb");
  if (!fp)
    return -1L;
  // Read problem dimensions: L, a, G
  if (fread(L, sizeof(unsigned), 1u, fp) != 1u)
    return -2L;
#ifndef NDEBUG
  (void)fprintf(stdout, "L = %u\n", *L);
  (void)fflush(stdout);
#endif /* !NDEBUG */
  if (fread(a, sizeof(unsigned), 1u, fp) != 1u)
    return -3L;
#ifndef NDEBUG
  (void)fprintf(stdout, "a = %u\n", *a);
  (void)fflush(stdout);
#endif /* !NDEBUG */
  if (fread(G, sizeof(unsigned), 1u, fp) != 1u)
    return -4L;
#ifndef NDEBUG
  (void)fprintf(stdout, "G = %u\n", *G);
  (void)fflush(stdout);
#endif /* !NDEBUG */
  // Read A and B buffers, *L * *a * *G * sizeof(double complex) each
  size_t size = (size_t)(*L) * (size_t)(*a) * (size_t)(*G);
  *A = (double complex*)calloc(size, sizeof(double complex));
  if (fread(*A, sizeof(double complex), size, fp) != size)
    return -5L;
  *B = (double complex*)calloc(size, sizeof(double complex));
  if (fread(*B, sizeof(double complex), size, fp) != size)
    return -6L;
  size = *a * (size_t)3u;
  *T_sizes = (unsigned*)calloc(size, sizeof(unsigned));
  *T = (double complex**)calloc(size, sizeof(double complex*));
  for (unsigned atom = 0u; atom < *a; ++atom) {
#ifndef NDEBUG
    (void)fprintf(stdout, "T[%u] = [ ", atom);
    (void)fflush(stdout);
#endif /* !NDEBUG */
    unsigned aa_size, ab_size, bb_size;
    // Taa
    if (fread(&aa_size, sizeof(unsigned), 1u, fp) != 1u)
      return (*a * -7L) - atom;
#ifndef NDEBUG
    (void)fprintf(stdout, "%u, ", aa_size);
    (void)fflush(stdout);
#endif /* !NDEBUG */
    (*T_sizes)[3u * atom + 0u] = aa_size;
    size = (size_t)aa_size * (size_t)aa_size;
    (*T)[3u * atom + 0u] = (double complex*)calloc(size, sizeof(double complex));
    double complex *const Taa = (*T)[3u * atom + 0u];
    if (fread(Taa, sizeof(double complex), size, fp) != size)
      return (*a * -8L) - atom;
    for (unsigned i = 0u; i < aa_size; ++i) {
      Taa[i * (size_t)aa_size + i] = creal(Taa[i * (size_t)aa_size + i]) + 0.0*I;
      for (unsigned j = i + 1u; j < aa_size; ++j)
        Taa[j * (size_t)aa_size + i] = conj(Taa[i * (size_t)aa_size + j]);
    }
    // Tab
    if (fread(&ab_size, sizeof(unsigned), 1u, fp) != 1u)
      return (*a * -9L) - atom;
#ifndef NDEBUG
    (void)fprintf(stdout, "%u, ", ab_size);
    (void)fflush(stdout);
#endif /* !NDEBUG */
    (*T_sizes)[3u * atom + 1u] = ab_size;
    size = (size_t)ab_size * (size_t)ab_size;
    (*T)[3u * atom + 1u] = (double complex*)calloc(size, sizeof(double complex));
    double complex *const Tab = (*T)[3u * atom + 1u];
    if (fread(Tab, sizeof(double complex), size, fp) != size)
      return (*a * -10L) - atom;
    // Tbb
    if (fread(&bb_size, sizeof(unsigned), 1u, fp) != 1u)
      return (*a * -11L) - atom;
#ifndef NDEBUG
    (void)fprintf(stdout, "%u ", bb_size);
    (void)fflush(stdout);
#endif /* !NDEBUG */
    (*T_sizes)[3u * atom + 2u] = bb_size;
    size = (size_t)bb_size * (size_t)bb_size;
    (*T)[3u * atom + 2u] = (double complex*)calloc(size, sizeof(double complex));
    double complex *const Tbb = (*T)[3u * atom + 2u];
    if (fread(Tbb, sizeof(double complex), size, fp) != size)
      return (*a * -12L) - atom;
    for (unsigned i = 0u; i < bb_size; ++i) {
      Tab[i * (size_t)bb_size + i] = creal(Tbb[i * (size_t)bb_size + i]) + 0.0*I;
      for (unsigned j = i + 1u; j < bb_size; ++j)
        Tab[j * (size_t)bb_size + i] = conj(Tbb[i * (size_t)bb_size + j]);
    }
#ifndef NDEBUG
    (void)fprintf(stdout, "]\n");
    (void)fflush(stdout);
#endif /* !NDEBUG */
  }
  // Buffers to compute norms
  size = (size_t)*a;
  *lmaxs = (unsigned*)calloc(size, sizeof(unsigned));
  if (fread(*lmaxs, sizeof(unsigned), size, fp) != size)
    return -13L;
#ifndef NDEBUG
  for (unsigned i = 0u; i < *a; ++i)
    (void)fprintf(stdout, "lmaxs[%u] = %u\n", i, (*lmaxs)[i]);
  (void)fflush(stdout);
#endif /* !NDEBUG */
  if (fread(max_lmax, sizeof(unsigned), 1u, fp) != 1u)
    return -14L;
#ifndef NDEBUG
  (void)fprintf(stdout, "max_lmax = %u\n", *max_lmax);
  (void)fflush(stdout);
#endif /* !NDEBUG */
  // u_norms: array of (*max_lmax+1) * *a doubles
  const unsigned max_lmax_1 = *max_lmax + 1u;
  size *= max_lmax_1;
  *u_norms = (double*)calloc(size, sizeof(double));
  if (fread(*u_norms, sizeof(double), size, fp) != size)
    return -15L;
#ifndef NDEBUG
  for (unsigned i = 0u; i < *a; ++i) {
    (void)fprintf(stdout, "u_norms[%u] = {", i);
    const size_t j = (size_t)i * max_lmax_1;
    const unsigned lmax = (*lmaxs)[i];
    for (unsigned l = 0u; l <= lmax; ++l) {
      const unsigned lm0 = l * l; // l^2
      const unsigned lme = lm0 + (l << 1u); // (l+1)^2 - 1
      (void)fprintf(stdout, " %#.16e %u %u", (*u_norms)[j + l], lm0, lme);
    }
    (void)fprintf(stdout, " }\n");
  }
  (void)fflush(stdout);
#endif /* !NDEBUG */
  const long ret = ftell(fp);
  // Done
  if (fclose(fp))
    return -16L;
  return ret;
}
