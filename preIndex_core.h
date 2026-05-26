#ifndef PREINDEX_CORE_H
#define PREINDEX_CORE_H

#define fragLength 10   // 조각(Read) 길이
#define fragNum 10      // 조각 개수
#define K_MER 3         // 자를 k-mer의 길이
#define HASH_SIZE (1 << (K_MER * 2))

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
SeedCandidate extractMax(MaxHeap* heap);
void freeMaxHeap(MaxHeap* heap);

void hashToKmer(int hash, char* out);
void printCountingIndex(const CountingIndex* index, char** frags);
void printTopSeeds(MaxHeap* heap, int topN);

#endif