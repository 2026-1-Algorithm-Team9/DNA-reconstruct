#ifndef PREINDEX_CORE_H
#define PREINDEX_CORE_H

#include <stddef.h>   // size_t

// ===== 데이터 생성 파라미터 =====
// refLength(원본 길이)만 바꾸면 리드 길이/개수가 커버리지에 맞춰 자동 조정됩니다.
// (권장: refLength >= 400. 너무 작으면 fragLength가 K_MER보다 작아져 동작 불가)
#define refLength    20000                                 // 원본 게놈 길이 (N)
#define COVERAGE     60                                    // 목표 커버리지 (배수)
#define fragLength   100                                   // 리드 길이 (L) — 현실적 짧은 리드
#define fragNum      ((refLength * COVERAGE) / fragLength) // 리드 개수 (M) = 커버리지 유지, 자동

// ===== 알고리즘 파라미터 =====
#define K_MER        10                                     // k-mer 길이 (시드 특이성)
#define HASH_SIZE    (1 << (K_MER * 2))                    // 해시 테이블 크기 (4^K)
#define MIN_OVERLAP  (fragLength / 2)                      // 조립 시 인정할 최소 겹침 길이
#define MAX_MISMATCH 2                                     // 겹침 검사 시 허용 미스매치 수

// ===== 시퀀싱 에러(미스매치) 주입 =====
// 리드를 떼올 때 각 염기를 이 확률(%)로 다른 염기로 치환(시퀀싱 에러/SNP 모사).
// 0이면 에러 없음. 발표용: 0 -> 1 -> 2 -> 5 로 바꿔가며 greedy vs consensus 비교.
#define ERROR_RATE_PERCENT 1

// ===== De Bruijn 그래프(consensus) 파라미터 =====
// 그래프용 k는 작아야 노드가 충분히 겹쳐 이어진다(보통 21 미만의 홀수). 여기선 13.
#define DBG_K        15                                    // De Bruijn k-mer 길이
#define DBG_HASH_SIZE (1 << (DBG_K * 2))                   // 4^DBG_K
// 빈도가 이 값 이하인 k-mer는 시퀀싱 에러로 간주해 그래프에서 제거.
// 커버리지가 높으면 올바른 k-mer는 빈도가 크고, 에러 k-mer는 1~2회뿐이다.
#define DBG_MIN_FREQ 2

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

// 5단계+: Consensus(다수결 투표) 보정 조립
// greedy 골격을 만든 뒤, 모든 리드를 골격에 매핑해 위치별로 A/C/G/T 투표,
// 다수결로 각 위치를 확정해 리드의 시퀀싱 에러를 커버리지로 눌러 정정한다.
char* assembleConsensus(const CountingIndex* index, const MaxHeap* heap, char** frags, int maxMismatch);

#endif
