#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

// ===============================================================
// 메모리 사용량 측정 (malloc 크기 합산, 간단 방식)
// ===============================================================
size_t countingIndexMemory(const CountingIndex* index) {
    if (index == NULL) return 0;

    size_t bytes = sizeof(CountingIndex);
    bytes += (size_t)HASH_SIZE * sizeof(int);                     // countArray
    bytes += (size_t)(HASH_SIZE + 1) * sizeof(int);              // prefixStart
    bytes += (size_t)HASH_SIZE * sizeof(int);                     // prefixEnd
    bytes += (size_t)index->totalKmers * sizeof(KmerOccurrence);  // occurrences
    return bytes;
}

size_t maxHeapMemory(const MaxHeap* heap) {
    if (heap == NULL) return 0;
    return sizeof(MaxHeap) + (size_t)heap->capacity * sizeof(SeedCandidate);
}

// ===============================================================
// 5단계: 우선순위 기반 조립 (Seed-and-Extend, de novo)
// ===============================================================

// 두 조각의 오버랩 구간을 비교하며 mismatch 수를 센다 (조기 종료)
int checkOverlapWithMismatch(const char* readSign, const char* targetRead, int overlapLen, int maxMismatch) {
    int mismatches = 0;
    for (int i = 0; i < overlapLen; i++) {
        if (readSign[i] != targetRead[i]) {
            mismatches++;
            if (mismatches > maxMismatch) return -1;   // 허용치 초과 시 즉시 실패
        }
    }
    return mismatches;
}

// MaxHeap의 고빈도 시드를 시작점으로, CountingIndex를 활용해 양방향으로 조각을 이어붙인다.
char* assembleReads(const CountingIndex* index, const MaxHeap* heap, char** frags, int maxMismatch) {
    if (heap->size == 0) return NULL;

    // 모든 조각이 일렬로 붙어도 넘치지 않을 만큼 큰 임시 버퍼 확보
    int maxPossibleLen = fragNum * fragLength * 2 + 1;
    char* tempBuffer = (char*)calloc(maxPossibleLen, sizeof(char));
    if (tempBuffer == NULL) return NULL;

    // 양방향 확장을 위해 시작 포인터를 버퍼 정중앙에 위치시킴
    int centerPos = maxPossibleLen / 2;
    char* assembledStart = tempBuffer + centerPos;

    // 가장 빈도 높은 1등 시드가 들어있는 첫 조각을 조립 시작점으로 선정
    int topHash = heap->data[0].hash;
    int startReadIdx = 0;
    if (index->prefixStart[topHash] < index->prefixEnd[topHash]) {
        startReadIdx = index->occurrences[index->prefixStart[topHash]].readIndex;
    }

    strcpy(assembledStart, frags[startReadIdx]);
    int assembledLen = (int)strlen(assembledStart);

    int* used = (int*)calloc(fragNum, sizeof(int));   // 중복 조립 방지용
    if (used != NULL) used[startReadIdx] = 1;

    // ---------- [방향 1] 오른쪽(뒤) 확장 ----------
    int progress = 1;
    while (progress) {
        progress = 0;
        int bestFragIdx = -1;
        int bestOverlap = 0;

        if (assembledLen >= K_MER) {
            // 현재 조립 서열의 끝 K_MER 글자를 해시로 변환
            char* tailKmerPtr = assembledStart + assembledLen - K_MER;
            int tailHash = 0;
            for (int j = 0; j < K_MER; j++)
                tailHash = (tailHash << 2) | charToBit(tailKmerPtr[j]);

            // 같은 끝 K_MER를 품은 후보 조각만 인덱스로 빠르게 추려서 검사
            int pStart = index->prefixStart[tailHash];
            int pEnd = index->prefixEnd[tailHash];
            for (int p = pStart; p < pEnd; p++) {
                int i = index->occurrences[p].readIndex;
                if (used && used[i]) continue;

                char* target = frags[i];
                // 오버랩 길이를 길게(fragLength)부터 줄여가며 가장 길게 겹치는 지점 탐색
                for (int len = fragLength; len >= MIN_OVERLAP; len--) {
                    if (assembledLen < len) continue;
                    char* assembledTail = assembledStart + assembledLen - len;
                    int mismatch = checkOverlapWithMismatch(assembledTail, target, len, maxMismatch);
                    if (mismatch != -1) {
                        if (len > bestOverlap) { bestOverlap = len; bestFragIdx = i; }
                        break;
                    }
                }
            }
        }

        if (bestFragIdx != -1 && bestOverlap > 0) {
            if (used) used[bestFragIdx] = 1;
            if ((assembledStart - tempBuffer) + assembledLen + (fragLength - bestOverlap) < maxPossibleLen) {
                strcat(assembledStart, frags[bestFragIdx] + bestOverlap);  // 겹친 뒤 꼬리만 결합
                assembledLen = (int)strlen(assembledStart);
                progress = 1;
            } else break;
        }
    }

    // ---------- [방향 2] 왼쪽(앞) 확장 ----------
    progress = 1;
    while (progress) {
        progress = 0;
        int bestFragIdx = -1;
        int bestOverlap = 0;

        if (assembledLen >= K_MER) {
            // 현재 조립 서열의 앞 K_MER 글자를 해시로 변환
            int headHash = 0;
            for (int j = 0; j < K_MER; j++)
                headHash = (headHash << 2) | charToBit(assembledStart[j]);

            int pStart = index->prefixStart[headHash];
            int pEnd = index->prefixEnd[headHash];
            for (int p = pStart; p < pEnd; p++) {
                int i = index->occurrences[p].readIndex;
                if (used && used[i]) continue;

                char* target = frags[i];
                for (int len = fragLength; len >= MIN_OVERLAP; len--) {
                    if (assembledLen < len) continue;
                    char* targetTail = target + fragLength - len;   // 후보 조각의 뒤쪽 len 글자
                    int mismatch = checkOverlapWithMismatch(targetTail, assembledStart, len, maxMismatch);
                    if (mismatch != -1) {
                        if (len > bestOverlap) { bestOverlap = len; bestFragIdx = i; }
                        break;
                    }
                }
            }
        }

        if (bestFragIdx != -1 && bestOverlap > 0) {
            if (used) used[bestFragIdx] = 1;
            int addLen = fragLength - bestOverlap;
            if (assembledStart - tempBuffer >= addLen) {
                assembledStart -= addLen;                            // 시작 포인터를 왼쪽으로 이동
                strncpy(assembledStart, frags[bestFragIdx], addLen); // 새로 확보된 앞 공간 채우기
                assembledLen += addLen;
                progress = 1;
            } else break;
        }
    }

    // 완성된 조립 서열만 떼어내 새 메모리에 복사
    char* finalAssembled = (char*)calloc(assembledLen + 1, sizeof(char));
    if (finalAssembled != NULL) {
        strncpy(finalAssembled, assembledStart, assembledLen);
        finalAssembled[assembledLen] = '\0';
    }

    if (used) free(used);
    free(tempBuffer);
    return finalAssembled;
}

// ===============================================================
// 5단계+: De Bruijn 그래프 기반 Consensus 조립 (에러 내성)
// ---------------------------------------------------------------
// 단순 greedy는 "리드를 통째로 이어붙여" 에러가 있으면 시드/겹침이 깨져
// 조립이 끊긴다. 여기서는 리드를 작은 k-mer(DBG_K)로 잘게 쪼개 그래프로 본다.
//
//   원리(커버리지로 에러를 이긴다):
//     - 올바른 k-mer는 여러 리드에 반복 등장 → 빈도 높음.
//     - 에러난 k-mer는 한두 번만 등장 → 빈도 낮음(DBG_MIN_FREQ 이하).
//     => 저빈도 k-mer를 버리면 '에러가 자동으로 걸러진' k-mer 집합이 남는다.
//
//   조립:
//     - k-mer 두 개가 (k-1)글자만큼 겹치면 한 글자씩 이어진다(De Bruijn 간선).
//     - 고빈도 시작 k-mer에서 출발해, 다음 글자(A/C/G/T) 중 빈도가 가장 높은
//       k-mer로 한 칸씩 전진하며 서열을 늘린다(= 위치별 다수결 consensus).
//     - 양방향(오른쪽/왼쪽)으로 확장해 원본을 복원한다.
// ===============================================================

// DBG_K 길이 k-mer 빈도표 생성 (리드를 슬라이딩하며 카운트). 에러 필터의 기반.
static int* buildDbgFreq(char** frags) {
    int* freq = (int*)calloc((size_t)DBG_HASH_SIZE, sizeof(int));
    if (freq == NULL) return NULL;
    int cut = (1 << (DBG_K * 2)) - 1;

    for (int i = 0; i < fragNum; i++) {
        char* read = frags[i];
        int h = 0;
        for (int j = 0; j < DBG_K; j++)
            h = (h << 2) | charToBit(read[j]);
        freq[h]++;
        for (int j = DBG_K; j < fragLength; j++) {
            h = ((h << 2) | charToBit(read[j])) & cut;
            freq[h]++;
        }
    }
    return freq;
}

// 현재 k-mer(curHash)의 뒤에 A/C/G/T를 붙여 만든 4개의 다음 k-mer 중
// 빈도가 가장 높은(>=DBG_MIN_FREQ) 것을 선택. 반환: 추가된 염기 비트(0~3), 없으면 -1.
static int bestNext(const int* freq, int curHash, int cut) {
    int bestBase = -1, bestFreq = DBG_MIN_FREQ;   // 임계값 이하(에러)는 무시
    for (int b = 0; b < 4; b++) {
        int nextHash = ((curHash << 2) | b) & cut;
        if (freq[nextHash] > bestFreq) { bestFreq = freq[nextHash]; bestBase = b; }
    }
    return bestBase;
}

// 현재 k-mer(curHash)의 앞에 A/C/G/T를 붙여 만든 4개의 이전 k-mer 중 최고 빈도 선택.
static int bestPrev(const int* freq, int curHash, int highBitShift) {
    int bestBase = -1, bestFreq = DBG_MIN_FREQ;
    for (int b = 0; b < 4; b++) {
        int prevHash = (curHash >> 2) | (b << highBitShift);
        if (freq[prevHash] > bestFreq) { bestFreq = freq[prevHash]; bestBase = b; }
    }
    return bestBase;
}

char* assembleConsensus(const CountingIndex* index, const MaxHeap* heap, char** frags, int maxMismatch) {
    (void)index; (void)heap; (void)maxMismatch;   // 이 방식은 DBG_K 빈도표만 사용

    int* freq = buildDbgFreq(frags);
    if (freq == NULL) return NULL;
    int cut = (1 << (DBG_K * 2)) - 1;
    int highBitShift = (DBG_K - 1) * 2;            // 앞쪽 염기를 끼울 비트 위치

    // 시작점: 빈도가 가장 높은 k-mer (가장 신뢰도 높은 지점에서 조립 시작)
    int seedHash = 0, seedFreq = -1;
    for (int h = 0; h < DBG_HASH_SIZE; h++) {
        if (freq[h] > seedFreq) { seedFreq = freq[h]; seedHash = h; }
    }
    if (seedFreq < DBG_MIN_FREQ) { free(freq); return NULL; }   // 쓸만한 k-mer 없음

    // 결과 버퍼: 양방향 확장을 위해 정중앙에서 시작
    int maxLen = refLength * 2 + fragLength + 16;   // 원본보다 넉넉히
    char* buf = (char*)calloc((size_t)maxLen, 1);
    if (buf == NULL) { free(freq); return NULL; }
    int center = maxLen / 2;
    char* start = buf + center;

    // 방문한 k-mer 재방문 방지(반복서열에서 무한루프 차단)
    char* visited = (char*)calloc((size_t)DBG_HASH_SIZE, 1);

    // 시드 k-mer를 가운데에 펼쳐 적기
    char seedStr[DBG_K + 1];
    { int h = seedHash; for (int i = DBG_K - 1; i >= 0; i--) { seedStr[i] = "ACGT"[h & 3]; h >>= 2; } seedStr[DBG_K] = '\0'; }
    memcpy(start, seedStr, DBG_K);
    int len = DBG_K;
    if (visited) visited[seedHash] = 1;

    // ---------- 오른쪽 확장 ----------
    int curHash = seedHash;
    while ((start - buf) + len < maxLen - 1) {
        int b = bestNext(freq, curHash, cut);
        if (b < 0) break;
        int nextHash = ((curHash << 2) | b) & cut;
        if (visited && visited[nextHash]) break;     // 반복 진입 차단
        start[len++] = "ACGT"[b];
        curHash = nextHash;
        if (visited) visited[nextHash] = 1;
    }

    // ---------- 왼쪽 확장 ----------
    curHash = seedHash;
    while (start - buf > 1) {
        int b = bestPrev(freq, curHash, highBitShift);
        if (b < 0) break;
        int prevHash = (curHash >> 2) | (b << highBitShift);
        if (visited && visited[prevHash]) break;
        start--; len++;
        start[0] = "ACGT"[b];
        curHash = prevHash;
        if (visited) visited[prevHash] = 1;
    }

    char* result = (char*)malloc((size_t)len + 1);
    if (result) { memcpy(result, start, len); result[len] = '\0'; }

    if (visited) free(visited);
    free(buf);
    free(freq);
    return result;
}


// ===============================================================
// 성능 리포트: 정확도(최적 오프셋 정렬) + 속도 + 메모리
// ===============================================================
void printPerformanceReport(const char* originalRef, const char* assembledRef, double duration, size_t memoryBytes) {
    if (originalRef == NULL || assembledRef == NULL) return;

    int N = (int)strlen(originalRef);
    int A = (int)strlen(assembledRef);

    // [정확도 측정] de novo 조립은 시작 위치/길이가 원본과 어긋나기 때문에
    // 0번부터 단순 비교하면 실력을 과소평가한다. 그래서 assembled를 original 위에서
    // 한 칸씩 밀어가며(offset) 일치 글자 수가 최대가 되는 위치를 찾아 그 값으로 채점한다.
    int bestMatches = 0;
    int bestOffset = 0;
    for (int d = -(A - 1); d <= N - 1; d++) {
        int matches = 0;
        for (int i = 0; i < A; i++) {
            int j = i + d;
            if (j >= 0 && j < N && assembledRef[i] == originalRef[j]) matches++;
        }
        if (matches > bestMatches) { bestMatches = matches; bestOffset = d; }
    }

    double accuracy = (N > 0) ? (100.0 * bestMatches / N) : 0.0;

    printf("\n================ [성능 분석 리포트] ================\n");
    printf("[정확도] 원본 복원율 : %.2f %% (%d / %d bp, 최적 오프셋 %d)\n",
           accuracy, bestMatches, N, bestOffset);
    printf("[속도]   알고리즘 시간: %.6f 초\n", duration);
    printf("[메모리] 사용량      : %zu bytes (%.2f KB)\n",
           memoryBytes, memoryBytes / 1024.0);
    printf("----------------------------------------------------\n");
    printf("원본 길이: %d | 조립 길이: %d\n", N, A);
    if (N <= 200) {
        printf("원본: %s\n", originalRef);
        printf("조립: %s\n", assembledRef);
    }
    printf("====================================================\n");
}