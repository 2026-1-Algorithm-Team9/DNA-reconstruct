#ifndef PREINDEX_CORE_H
#define PREINDEX_CORE_H

#include <stddef.h>   // size_t

// ===== 데이터 생성 파라미터 =====
// refLength(원본 길이)만 바꾸면 리드 길이/개수가 커버리지에 맞춰 자동 조정됩니다.
// (권장: refLength >= 400. 너무 작으면 fragLength가 K_MER보다 작아져 동작 불가)
#define refLength    100000                                  // 원본 게놈 길이 (N)
#define COVERAGE     30                                    // 목표 커버리지 (배수)
#define fragLength   (refLength / 20)                      // 리드 길이 (L) = 원본의 1/20, 자동
#define fragNum      ((refLength * COVERAGE) / fragLength) // 리드 개수 (M) = 커버리지 유지, 자동

// ===== 알고리즘 파라미터 =====
#define K_MER        10                                     // k-mer 길이 (시드 특이성)
#define HASH_SIZE    (1 << (K_MER * 2))                    // 해시 테이블 크기 (4^K)
#define MIN_OVERLAP  (fragLength / 2)                      // 조립 시 인정할 최소 겹침 길이
#define MAX_MISMATCH 2                                     // 겹침 검사 시 허용 미스매치 수

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

// 1~3단계: Counting Sort 기반 압축 인덱스
CountingIndex* buildCountingIndex(char** frags);
void freeCountingIndex(CountingIndex* index);
size_t countingIndexMemory(const CountingIndex* index);

// 4단계: Max Heap 시드
MaxHeap* buildSeedHeap(const CountingIndex* index);
SeedCandidate extractMax(MaxHeap* heap);
void freeMaxHeap(MaxHeap* heap);
size_t maxHeapMemory(const MaxHeap* heap);

// 유틸 / 출력
void hashToKmer(int hash, char* out);
void printCountingIndex(const CountingIndex* index, char** frags);
void printTopSeeds(MaxHeap* heap, int topN);

// 5단계: 우선순위 기반 조립 + 성능 리포트
char* assembleReads(const CountingIndex* index, const MaxHeap* heap, char** frags, int maxMismatch);
void printPerformanceReport(const char* originalRef, const char* assembledRef, double duration, size_t memoryBytes);

#endif
