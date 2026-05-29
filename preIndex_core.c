#include <stdio.h>
#include <stdlib.h>

#include "preIndex_core.h"

static int charToBit(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
    }

    return 0;
}

static int* freqArray(char** frags) {
    int* countArray = (int*)calloc(HASH_SIZE, sizeof(int));
    if (countArray == NULL) return NULL;

    int cut = (1 << (K_MER * 2)) - 1;

    for (int i = 0; i < fragNum; i++) {
        char* read = frags[i];
        int currentHash = 0;

        for (int j = 0; j < K_MER; j++) {
            currentHash = (currentHash << 2) | charToBit(read[j]);
        }
        countArray[currentHash]++;

        for (int j = K_MER; j < fragLength; j++) {
            currentHash = ((currentHash << 2) | charToBit(read[j])) & cut;
            countArray[currentHash]++;
        }
    }

    return countArray;
}

static char bitToChar(int value) {
    switch (value) {
        case 0: return 'A';
        case 1: return 'C';
        case 2: return 'G';
        case 3: return 'T';
    }

    return 'A';
}

void hashToKmer(int hash, char* out) {
    for (int i = K_MER - 1; i >= 0; i--) {
        out[i] = bitToChar(hash & 3);
        hash >>= 2;
    }
    out[K_MER] = '\0';
}

CountingIndex* buildCountingIndex(char** frags) {
    CountingIndex* index = (CountingIndex*)malloc(sizeof(CountingIndex));
    if (index == NULL) return NULL;

    index->countArray = freqArray(frags);
    if (index->countArray == NULL) {
        free(index);
        return NULL;
    }

    index->totalKmers = fragNum * (fragLength - K_MER + 1);
    index->prefixStart = (int*)malloc((HASH_SIZE + 1) * sizeof(int));
    index->prefixEnd = (int*)malloc(HASH_SIZE * sizeof(int));
    index->occurrences = (KmerOccurrence*)malloc(index->totalKmers * sizeof(KmerOccurrence));

    if (index->prefixStart == NULL || index->prefixEnd == NULL || index->occurrences == NULL) {
        free(index->occurrences);
        free(index->prefixEnd);
        free(index->prefixStart);
        free(index->countArray);
        free(index);
        return NULL;
    }

    index->prefixStart[0] = 0;
    for (int hash = 1; hash <= HASH_SIZE; hash++) {
        index->prefixStart[hash] = index->prefixStart[hash - 1] + index->countArray[hash - 1];
    }

    for (int hash = 0; hash < HASH_SIZE; hash++) {
        index->prefixEnd[hash] = index->prefixStart[hash] + index->countArray[hash];
    }

    int* currentOffset = (int*)calloc(HASH_SIZE, sizeof(int));
    if (currentOffset == NULL) {
        free(index->occurrences);
        free(index->prefixEnd);
        free(index->prefixStart);
        free(index->countArray);
        free(index);
        return NULL;
    }

    int cut = (1 << (K_MER * 2)) - 1;

    for (int readIndex = 0; readIndex < fragNum; readIndex++) {
        char* read = frags[readIndex];
        int currentHash = 0;

        for (int j = 0; j < K_MER; j++) {
            currentHash = (currentHash << 2) | charToBit(read[j]);
        }

        int firstPosition = index->prefixStart[currentHash] + currentOffset[currentHash];
        index->occurrences[firstPosition].readIndex = readIndex;
        index->occurrences[firstPosition].offset = 0;
        currentOffset[currentHash]++;

        for (int j = K_MER; j < fragLength; j++) {
            currentHash = ((currentHash << 2) | charToBit(read[j])) & cut;

            int offset = j - K_MER + 1;
            int position = index->prefixStart[currentHash] + currentOffset[currentHash];

            index->occurrences[position].readIndex = readIndex;
            index->occurrences[position].offset = offset;
            currentOffset[currentHash]++;
        }
    }

    free(currentOffset);
    return index;
}

void freeCountingIndex(CountingIndex* index) {
    if (index == NULL) return;

    free(index->occurrences);
    free(index->prefixEnd);
    free(index->prefixStart);
    free(index->countArray);
    free(index);
}

static void swapSeed(SeedCandidate* a, SeedCandidate* b) {
    SeedCandidate temp = *a;
    *a = *b;
    *b = temp;
}

static void heapifyUp(MaxHeap* heap, int index) {
    while (index > 0) {
        int parent = (index - 1) / 2;

        if (heap->data[parent].frequency >= heap->data[index].frequency) break;

        swapSeed(&heap->data[parent], &heap->data[index]);
        index = parent;
    }
}

static void heapifyDown(MaxHeap* heap, int index) {
    while (1) {
        int left = index * 2 + 1;
        int right = index * 2 + 2;
        int largest = index;

        if (left < heap->size && heap->data[left].frequency > heap->data[largest].frequency) {
            largest = left;
        }

        if (right < heap->size && heap->data[right].frequency > heap->data[largest].frequency) {
            largest = right;
        }

        if (largest == index) break;

        swapSeed(&heap->data[index], &heap->data[largest]);
        index = largest;
    }
}

MaxHeap* buildSeedHeap(const CountingIndex* index) {
    MaxHeap* heap = (MaxHeap*)malloc(sizeof(MaxHeap));
    if (heap == NULL) return NULL;

    heap->data = (SeedCandidate*)malloc(HASH_SIZE * sizeof(SeedCandidate));
    heap->size = 0;
    heap->capacity = HASH_SIZE;

    if (heap->data == NULL) {
        free(heap);
        return NULL;
    }

    for (int hash = 0; hash < HASH_SIZE; hash++) {
        if (index->countArray[hash] == 0) continue;

        heap->data[heap->size].hash = hash;
        heap->data[heap->size].frequency = index->countArray[hash];
        heapifyUp(heap, heap->size);
        heap->size++;
    }

    return heap;
}

SeedCandidate extractMax(MaxHeap* heap) {
    SeedCandidate empty = {-1, 0};
    if (heap == NULL || heap->size == 0) return empty;

    SeedCandidate top = heap->data[0];
    heap->size--;
    heap->data[0] = heap->data[heap->size];
    heapifyDown(heap, 0);

    return top;
}

void freeMaxHeap(MaxHeap* heap) {
    if (heap == NULL) return;

    free(heap->data);
    free(heap);
}

void printCountingIndex(const CountingIndex* index, char** frags) {
    char kmer[K_MER + 1];

    printf("=== 단계 3: Counting Sort 기반 압축 인덱스 ===\n");
    printf("총 k-mer 개수: %d\n", index->totalKmers);

    for (int hash = 0; hash < HASH_SIZE; hash++) {
        if (index->countArray[hash] == 0) continue;

        hashToKmer(hash, kmer);
        printf("해시 [%2d] k-mer [%s] -> 시작 %2d, 끝 %2d, 빈도 %d\n",
               hash,
               kmer,
               index->prefixStart[hash],
               index->prefixEnd[hash] - 1,
               index->countArray[hash]);

        for (int pos = index->prefixStart[hash]; pos < index->prefixEnd[hash]; pos++) {
            KmerOccurrence occurrence = index->occurrences[pos];
            printf("    압축배열[%2d] = read %d, offset %d, 조각 내 k-mer %.*s\n",
                   pos,
                   occurrence.readIndex,
                   occurrence.offset,
                   K_MER,
                   frags[occurrence.readIndex] + occurrence.offset);
        }
    }

    printf("\n");
}

void printTopSeeds(MaxHeap* heap, int topN) {
    char kmer[K_MER + 1];

    printf("=== 단계 4: Max Heap 기반 상위 시드 추출 ===\n");
    for (int rank = 1; rank <= topN && heap->size > 0; rank++) {
        SeedCandidate candidate = extractMax(heap);
        hashToKmer(candidate.hash, kmer);
        printf("Top %d Seed -> hash %2d, k-mer [%s], 빈도 %d\n",
               rank,
               candidate.hash,
               kmer,
               candidate.frequency);
    }
    printf("\n");
}