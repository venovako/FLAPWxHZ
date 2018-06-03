static long split_data
(const char *const file,
 const unsigned L,
 const unsigned a,
 const unsigned G,
 const unsigned max_lmax,
 const unsigned *const T_sizes,
 const double complex *const A,
 const double complex *const B,
 const double complex *const *const T,
 const unsigned *const lmaxs,
 const double *const u_norms)
{
  const size_t fnl = strlen(file) + 3u;
  char fn[fnl];
  (void)snprintf(fn, fnl, "%s.X", file);
  FILE *const fX = fopen(fn, "wb");
  if (!fX)
    return -1L;
  (void)snprintf(fn, fnl, "%s.T", file);
  FILE *const fT = fopen(fn, "wb");
  if (!fT)
    return -2L;
  (void)snprintf(fn, fnl, "%s.U", file);
  FILE *const fU = fopen(fn, "wb");
  if (!fU)
    return -3L;
  double complex *const TF = (double complex*)calloc(4UL * L * L, sizeof(double complex));
  if (!TF)
    return -4L;
  double *const U = (double*)calloc(L, sizeof(double));
  if (!U)
    return -5L;
  const size_t ML = max_lmax + 1UL;
  for (unsigned atom = 0u; atom < a; ++atom) {
    for (unsigned j = 0u; j < G; ++j) {
      if (fwrite(A + ((size_t)j * a + atom) * L, sizeof(double complex), L, fX) != L)
        return -6L;
      if (fwrite(B + ((size_t)j * a + atom) * L, sizeof(double complex), L, fX) != L)
        return -7L;
    }
    const size_t aa_size = T_sizes[3u * atom + 0u];
    if (aa_size != L)
      return -8L;
    const size_t bb_size = T_sizes[3u * atom + 2u];
    if (bb_size != L)
      return -9L;
    const size_t ab_size = T_sizes[3u * atom + 1u];
    if (ab_size != L)
      return -10L;
    const size_t T_size1 = aa_size + bb_size;
    const size_t T_size2 = T_size1 * T_size1;
    const double complex *const Taa = T[3u * atom + 0u];
    // tril(Taa) -> TF
    for (unsigned j = 0u; j < (unsigned)aa_size; ++j) {
      TF[j * T_size1 + j] = creal(Taa[j * aa_size + j]);
      for (unsigned i = j + 1u; i < (unsigned)aa_size; ++i)
        TF[j * T_size1 + i] = Taa[j * aa_size + i];
    }
    for (unsigned j = 0u; j < (unsigned)aa_size; ++j)
      for (unsigned i = 0u; i < j; ++i)
        TF[j * T_size1 + i] = conj(TF[i * T_size1 + j]);
    const double complex *const Tbb = T[3u * atom + 2u];
    // tril(Tbb) -> TF
    for (unsigned j = 0u; j < (unsigned)bb_size; ++j) {
      TF[(j + aa_size) * T_size1 + (j + aa_size)] = creal(Tbb[j * bb_size + j]);
      for (unsigned i = j + 1u; i < (unsigned)bb_size; ++i)
        TF[(j + aa_size) * T_size1 + (i + aa_size)] = Tbb[j * bb_size + i];
    }
    for (unsigned j = aa_size; j < (unsigned)T_size1; ++j)
      for (unsigned i = aa_size; i < j; ++i)
        TF[j * T_size1 + i] = conj(TF[i * T_size1 + j]);
    const double complex *const Tab = T[3u * atom + 1u];
    // Tab -> TF
    for (unsigned j = 0u; j < (unsigned)bb_size; ++j)
      for (unsigned i = 0u; i < (unsigned)aa_size; ++i)
        TF[(j + aa_size) * T_size1 + i] = Tab[j * aa_size + i];
    // Tab^H -> TF
    for (unsigned j = aa_size; j < (unsigned)T_size1; ++j)
      for (unsigned i = 0u; i < (unsigned)aa_size; ++i)
        TF[i * T_size1 + j] = conj(TF[j * T_size1 + i]);
    if (fwrite(TF, sizeof(double complex), T_size2, fT) != T_size2)
      return -11L;
    const size_t atML = atom * ML;
    const unsigned lmax = lmaxs[atom];
    for (unsigned l = 0u; l <= lmax; ++l) {
      const double u = u_norms[atML + l];
      const unsigned lm0 = l * l; // l^2
      const unsigned lme = lm0 + (l << 1u); // (l+1)^2 - 1
      for (unsigned lm = lm0; lm <= lme; ++lm)
        U[lm] = u;
    }
    if (fwrite(U, sizeof(double), L, fU) != L)
      return -12L;
  }
  if (fclose(fU))
    return -13L;
  if (fclose(fT))
    return -14L;
  if (fclose(fX))
    return -15L;
  free(U);
  free(TF);
  (void)fprintf(stdout, "L = %u\na = %u\nG = %u\n", L, a, G);
  (void)fflush(stdout);
  return (long)a;
}
