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

void printTopSeeds(const MaxHeap* heap, int topN) {
    printf("=== 상위 %d개 고빈도 Seed (K-mer) ===\n", topN);
    int limit = (heap->size < topN) ? heap->size : topN;
    for (int i = 0; i < limit; i++) {
        char kmerStr[K_MER + 1];
        hashToKmer(heap->data[i].hash, kmerStr);
        printf("[%d등] K-mer: %s | 등장 빈도: %d회\n", i + 1, kmerStr, heap->data[i].frequency);
    }
    printf("\n");
}

// ===============================================================
// 조립 테스트용 코드
// ===============================================================

/*
 * 1. 두 조각의 오버랩 구간을 비교하며 mismatch 개수를 세는 함수 (조기 종료 적용)
 */
int checkOverlapWithMismatch(const char* readSign, const char* targetRead, int overlapLen, int maxMismatch) {
    int mismatches = 0;
    
    for (int i = 0; i < overlapLen; i++) {
        if (readSign[i] != targetRead[i]) {
            mismatches++;
            
            // [조기 종료] 허용치를 넘는 순간 다음 글자는 보지도 않고 즉시 실패 리턴
            if (mismatches > maxMismatch) {
                return -1; 
            }
        }
    }
    return mismatches; // 매칭 성공 시 실제 발생한 mismatch 수 반환
}

/*
 * 2. Max Heap 시드와 CountingIndex를 활용해 조각들을 초고속 조립하는 메인 함수 (방법 A 적용)
 */
/*
 * CountingIndex(해시 테이블)와 MaxHeap(고빈도 시드)을 활용해 
 * 염기서열 조각들을 초고속 및 99% 고정밀로 양방향 조립하는 메인 함수
 */
char* assembleReads(const CountingIndex* index, const MaxHeap* heap, char** frags, int maxMismatch) {
    // 힙에 데이터(시드)가 없으면 조립을 시작할 수 없으므로 즉시 종료
    if (heap->size == 0) return NULL;

    // 1. 모든 조각이 겹치지 않고 일렬로 다 붙어도 차고 넘치도록 거대한 임시 메모리 버퍼 생성
    int maxPossibleLen = fragNum * fragLength * 2 + 1;
    char* tempBuffer = (char*)calloc(maxPossibleLen, sizeof(char));
    if (tempBuffer == NULL) return NULL;

    // 2. 양방향(왼쪽/오른쪽) 확장을 자유롭게 하기 위해 포인터를 거대 버퍼의 '정중앙'에 위치시킴
    int centerPos = maxPossibleLen / 2;
    char* assembledStart = tempBuffer + centerPos;

    // 3. MaxHeap에서 등장 빈도가 가장 높았던 1등 대장 해시(시드)를 추출
    int topHash = heap->data[0].hash;
    int startReadIdx = -1;
    int startPos = index->prefixStart[topHash];
    int endPos = index->prefixEnd[topHash];
    
    // 1등 시드가 포함된 첫 번째 기준 조각(Read 번호)을 인덱스 지도에서 조회
    if (startPos != -1 && startPos < endPos) {
        startReadIdx = index->occurrences[startPos].readIndex;
    } else {
        startReadIdx = 0; // 안전망: 예외 발생 시 0번째 조각을 기준점으로 삼음
    }

    // 선정된 첫 번째 조각을 정중앙 시작 서열로 복사
    strcpy(assembledStart, frags[startReadIdx]);
    int assembledLen = strlen(assembledStart);

    // 중복 조립을 방지하기 위한 '이미 사용된 조각 추적 배열' 생성 및 첫 조각 체크
    int* used = (int*)calloc(fragNum, sizeof(int));
    if (used != NULL) used[startReadIdx] = 1;

    // =================================================================
    // [방향 1] 오른쪽(뒤쪽) 확장 루프 (인덱스 기반 초고속 유동 탐색)
    // =================================================================
    int progress = 1;
    while (progress) {
        progress = 0; // 이번 회차에 새로 붙인 조각이 없다면 루프가 멈추도록 초기화
        int bestFragIdx = -1;  // 가장 잘 어울리는 후보 조각 번호 저장
        int bestOverlap = 0;   // 가장 길게 겹치는 오버랩 길이 저장

        // 현재까지 조립된 전체 서열의 길이가 최소 K_MER(6글자) 이상일 때만 해싱 탐색 가능
        if (assembledLen >= K_MER) {
            // 현재 조립 서열의 '맨 오른쪽 끝 6글자'를 가리키는 포인터 추출
            char* tailKmerPtr = assembledStart + assembledLen - K_MER;
            int tailHash = 0;
            
            // 끝자락 6글자 염기를 비트 연산을 통해 정수 해시 값(tailHash)으로 변환
            for (int j = 0; j < K_MER; j++) {
                tailHash = (tailHash << 2) | charToBit(tailKmerPtr[j]);
            }

            // 인덱스 맵에서 이 6글자 해시를 품고 있는 후보 조각들의 명단 구간(pStart ~ pEnd) 확보
            int pStart = index->prefixStart[tailHash];
            int pEnd = index->prefixEnd[tailHash];

            // [속도 혁신 핵심]: 모든 조각(210개)을 뒤지지 않고, 방금 솎아낸 소수의 후보 조각만 순회
            for (int p = pStart; p < pEnd; p++) {
                int i = index->occurrences[p].readIndex; // 후보 조각 번호 꺼내기
                if (used && used[i]) continue;           // 이미 조립에 써먹은 조각이면 패스

                char* target = frags[i];

                // [정확도 혁신 핵심]: 무작위 자르기로 인한 오프셋 어긋남을 해결하기 위해 
                // 후보 조각 내에서 오버랩 길이(len)를 50부터 최소 기준인 15까지 줄여가며 정밀 대조
                for (int len = fragLength; len >= MIN_OVERLAP; len--) {
                    if (assembledLen < len) continue; // 포갤 서열이 검사하려는 길이보다 짧으면 패스

                    // 기존 서열의 오른쪽 끝 len 글자와 후보 조각의 앞 len 글자 매칭 준비
                    char* assembledTail = assembledStart + assembledLen - len;
                    int mismatch = checkOverlapWithMismatch(assembledTail, target, len, maxMismatch);

                    // 허용된 미스매치 범위 내에서 매칭 성공 시
                    if (mismatch != -1) {
                        // 기존에 찾은 매칭보다 더 길게 겹치는 최적의 조각인 경우 정보 갱신
                        if (len > bestOverlap) {
                            bestOverlap = len;
                            bestFragIdx = i;
                        }
                        break; // 이 조각에서 가장 길게 포개지는 부위를 찾았으므로 len 감소 루프 탈출
                    }
                }
            }
        }

        // 인덱스 내부 정밀 탐색을 통해 우측에 이어 붙일 최고의 짝꿍 조각을 찾았다면 결합 진행
        if (bestFragIdx != -1 && bestOverlap > 0) {
            if (used) used[bestFragIdx] = 1; // 해당 조각 사용 완료 처리
            
            // 거대 버퍼의 남은 공간을 넘지 않는지 안전성 검사 후 이어 붙이기
            if ((assembledStart - tempBuffer) + assembledLen + (fragLength - bestOverlap) < maxPossibleLen) {
                // 겹치는 부위(bestOverlap) 이후의 꼬리 문자열만 기존 서열 맨 뒤에 결합(strcat)
                strcat(assembledStart, frags[bestFragIdx] + bestOverlap);
                assembledLen = strlen(assembledStart); // 조립 서열 길이 새로 갱신
                progress = 1;                          // 조각을 찾았으니 다음 회차 while 루프 연장
            } else {
                break;
            }
        }
    }

    // =================================================================
    // [방향 2] 왼쪽(앞쪽) 확장 루프 (인덱스 기반 초고속 유동 탐색)
    // =================================================================
    progress = 1;
    while (progress) {
        progress = 0; // 초기화
        int bestFragIdx = -1;
        int bestOverlap = 0;

        if (assembledLen >= K_MER) {
            // 이번에는 현재 조립 서열의 '맨 왼쪽 시작(머리맡) 6글자' 해시 값 구하기
            int headHash = 0;
            for (int j = 0; j < K_MER; j++) {
                headHash = (headHash << 2) | charToBit(assembledStart[j]);
            }

            // 인덱스 맵에서 맨 앞 6글자를 품고 있는 좌측 확장용 후보 조각 명단 확보
            int pStart = index->prefixStart[headHash];
            int pEnd = index->prefixEnd[headHash];

            // 솎아내진 좌측 확장 후보 조각들만 제한적으로 스캔
            for (int p = pStart; p < pEnd; p++) {
                int i = index->occurrences[p].readIndex;
                if (used && used[i]) continue;

                char* target = frags[i];

                // 후보 조각 안에서 왼쪽으로 포갤 길이(len)를 50부터 15까지 줄여가며 검사
                for (int len = fragLength; len >= MIN_OVERLAP; len--) {
                    if (assembledLen < len) continue;

                    // 후보 조각의 맨 뒤쪽 len 글자와 기존 조립 서열의 맨 앞 len 글자 매칭 준비
                    char* targetTail = target + fragLength - len;
                    int mismatch = checkOverlapWithMismatch(targetTail, assembledStart, len, maxMismatch);

                    if (mismatch != -1) {
                        if (len > bestOverlap) {
                            bestOverlap = len;
                            bestFragIdx = i;
                        }
                        break; // 가장 길게 포개지는 부위를 찾았으므로 탈출
                    }
                }
            }
        }

        // 왼쪽으로 이어 붙일 최고의 조각을 찾았다면 결합 진행
        if (bestFragIdx != -1 && bestOverlap > 0) {
            if (used) used[bestFragIdx] = 1;

            // 새로 왼쪽에 덧붙여질 늘어날 길이 계산 (50 - 겹친길이)
            int addLen = fragLength - bestOverlap;
            
            // 거대 버퍼의 왼쪽 빈 공간이 충분한지 확인 후 포인터 조작을 통한 결합
            if (assembledStart - tempBuffer >= addLen) {
                assembledStart -= addLen; // ★시작점 포인터를 왼쪽으로 addLen만큼 전진 이동★
                
                // 이동해서 확보된 왼쪽 새 공간에 후보 조각의 앞부분 내용을 복사(strncpy)하여 채워넣음
                strncpy(assembledStart, frags[bestFragIdx], addLen);
                assembledLen += addLen; // 조립 서열 길이 새로 증가
                progress = 1;           // 다음 루프 연장
            } else {
                break;
            }
        }
    }

    // 4. 거대 임시 버퍼에서 마침내 완성된 순수 조립 서열 문자열만 똑 떼어내어 새 공간에 할당
    char* finalAssembled = (char*)calloc(assembledLen + 1, sizeof(char));
    if (finalAssembled != NULL) {
        strncpy(finalAssembled, assembledStart, assembledLen); // 알맹이만 복사
        finalAssembled[assembledLen] = '\0'; // 문자열 종료 문자 삽입
    }

    // 5. 사용이 끝난 추적 배열 및 임시 거대 버퍼 메모리 반환 후 최종 서열 리턴
    if (used) free(used);
    free(tempBuffer);

    return finalAssembled; // 최종 결과물 주소 반환
}


/**
 * 5단계: 최종 조립 결과 및 정확도/속도 검증 출력 함수
 */
void printPerformanceReport(char* originalRef, char* assembledRef, double duration) {
    if (!originalRef || !assembledRef) return;

    int origLen = strlen(originalRef);
    int asmbLen = strlen(assembledRef);

    printf("\n================ [성능 분석 리포트] ================\n");
    printf("▶ 총 실행 시간  : %.6f 초\n", duration);
    printf("▶ 원본 서열 길이: %d\n", origLen);
    printf("▶ 조립 서열 길이: %d\n", asmbLen);
    printf("----------------------------------------------------\n");
    
    printf("원본: %s\n", originalRef);
    printf("조립: %s\n", assembledRef);
    printf("----------------------------------------------------\n");

    int matchCount = 0;
    double accuracy = 0.0;

    // strstr을 사용하여 원본 서열 안에 조립 서열이 통째로 예쁘게 들어가는지 검사
    if (strstr(originalRef, assembledRef) != NULL) {
        matchCount = asmbLen;
        accuracy = ((double)matchCount / origLen) * 100.0;
        printf("[성능 분석 리포트 - 부분 일치 성공]\n");
    } else {
        // 단방향 확장 특성상 시작점이 밀릴 수 있으므로, 원본에서 가장 잘 매칭되는 오프셋을 찾아서 비교하면 좋습니다.
        // 여기서는 기존 1대1 매칭 로직의 안전망을 유지합니다.
        int checkLen = (origLen < asmbLen) ? origLen : asmbLen;
        for (int i = 0; i < checkLen; i++) {
            if (originalRef[i] == assembledRef[i]) {
                matchCount++;
            }
        }
        accuracy = ((double)matchCount / origLen) * 100.0;
        printf("[성능 분석 리포트 - 위치 어긋남 또는 Mismatch 포함 복원]\n");
    }

    printf("▶ 최종 매칭 염기 수: %d / %d\n", matchCount, origLen);
    printf("▶ 원본 서열 복원율  : %.2f %%\n", accuracy);
    printf("====================================================\n");
}