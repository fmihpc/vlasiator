#ifdef __GNUC__

inline void prefetch(const char *addr, int sz){
   const int Nlines = sz/CACHELINE_BYTES+1;
   for (int c=0; c<Nlines; c++)
      __builtin_prefetch(addr+CACHELINE_BYTES*c);
}

#else

inline void prefetch(const char*,int);

#endif



prefetch(reinterpret_cast<char*>(&listptr[1]),sizeof(InternalParticleBlock));




