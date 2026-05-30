#ifndef PREINDEX_CORE_H
#define PREINDEX_CORE_H

#define fragLength 50   // 조각(Read) 길이
#define fragNum 210     // 조각 개수
#define K_MER 6         // 자를 k-mer의 길이
#define HASH_SIZE (1 << (K_MER * 2))

#define MAX_MISMATCH 3       // 허용할 최대 mismatch 수
#define MIN_OVERLAP 15        // 조각 간에 최소한으로 겹쳐야 하는 글자 수

typedef struct {
    int readIndex;
    int offset;
} KmerOccurrence;

typedef struct {
    int* countArray;
    int* prefixStart;
    int* prefixEnd;
    KmerOccurrence* occurrences;
    int totalKmers;
} CountingIndex;

typedef struct {
    int hash;
    int frequency;
} SeedCandidate;

typedef struct {
    SeedCandidate* data;
    int size;
    int capacity;
} MaxHeap;

CountingIndex* buildCountingIndex(char** frags);
void freeCountingIndex(CountingIndex* index);

MaxHeap* buildSeedHeap(const CountingIndex* index);
void freeMaxHeap(MaxHeap* heap);

void hashToKmer(int hash, char* out);
void printCountingIndex(const CountingIndex* index, char** frags);
void printTopSeeds(const MaxHeap* heap, int topN);

int checkOverlapWithMismatch(const char* readSign, const char* targetRead, int overlapLen, int maxMismatch);
char* assembleReads(const CountingIndex* index, const MaxHeap* heap, char** frags, int maxMismatch);
void printPerformanceReport(char* originalRef, char* assembledRef, double duration);

#endif