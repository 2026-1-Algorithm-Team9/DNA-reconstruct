#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "preIndex_core.h"

// OS와 무관하게 충분히 큰 난수 생성.
// Windows(MinGW)의 RAND_MAX는 32767로 작아, rand()만 쓰면 startIndex가
// 원본 앞부분에만 몰려 뒤쪽이 리드로 안 덮인다(복원율 급락). rand()를 두 번
// 조합해 약 30비트 난수를 만들어 어느 환경에서도 전 범위가 고르게 나오게 한다.
static long bigRand(void) {
    return ((long)rand() << 15) | (long)rand();
}

char* makeRef() {
    char basis[4] = {'A', 'T', 'C', 'G'};

    // 원본 염기서열은 A/C/G/T를 완전 랜덤으로 생성 (가상의 정답 게놈)
    char* ref = (char*)malloc((refLength + 1) * sizeof(char));
    if (ref == NULL) return NULL;
    srand((unsigned int)time(NULL));

    for (int i = 0; i < refLength; i++) {
        int index = rand() % 4;
        ref[i] = basis[index];
    }
    ref[refLength] = '\0';

    return ref;
}

// 한 염기를 자신과 다른 3개 중 하나로 무작위 치환 (= 미스매치 1개 발생)
static char mutateBase(char original) {
    char basis[4] = {'A', 'C', 'G', 'T'};
    char c;
    do {
        c = basis[rand() % 4];
    } while (c == original);
    return c;
}

char** makeFrag(char* ref, int* outInjectedErrors) {
    char** frags = (char**)malloc(fragNum * sizeof(char*));
    if (frags == NULL) return NULL;

    int injected = 0;

    for (int i = 0; i < fragNum; i++) {
        frags[i] = (char*)malloc((fragLength + 1) * sizeof(char));

        int startIndex = bigRand() % (refLength - fragLength + 1);   // 원본에서 랜덤 시작 위치 (OS 무관)
        for (int j = 0; j < fragLength; j++) {
            char base = ref[startIndex + j];                         // 원본의 substring을 떼옴

            // ERROR_RATE_PERCENT 확률로 시퀀싱 에러(미스매치) 주입
            if (ERROR_RATE_PERCENT > 0 && (rand() % 100) < ERROR_RATE_PERCENT) {
                base = mutateBase(base);
                injected++;
            }
            frags[i][j] = base;
        }
        frags[i][fragLength] = '\0';
    }

    if (outInjectedErrors != NULL) *outInjectedErrors = injected;
    return frags;
}

// 메모리 사용량 합산 (malloc 기준 간단 측정)
static size_t calcMemory(const CountingIndex* idx, const MaxHeap* heap, size_t assembledLen) {
    size_t bytes = 0;
    bytes += countingIndexMemory(idx);
    bytes += maxHeapMemory(heap);
    bytes += (size_t)fragNum * sizeof(char*) + (size_t)fragNum * (fragLength + 1);  // 리드
    bytes += (size_t)(refLength + 1);                                               // 원본
    bytes += (size_t)(fragNum * fragLength * 2 + 1);                                // 조립 임시버퍼(peak)
    bytes += assembledLen + 1;                                                      // 결과 서열
    return bytes;
}

int main(void) {
    int injectedErrors = 0;
    char* ref = makeRef();
    char** frags = (ref != NULL) ? makeFrag(ref, &injectedErrors) : NULL;

    if (ref == NULL || frags == NULL) {
        printf("초기 데이터 생성에 실패했습니다.\n");
        free(frags);
        free(ref);
        return 1;
    }

    long totalBases = (long)fragNum * fragLength;
    printf("=== 데이터 생성 파라미터 ===\n");
    printf("원본 길이(N): %d | 리드 길이(L): %d | 리드 개수(M): %d | 커버리지: %.1f배 | K-mer: %d\n",
           refLength, fragLength, fragNum,
           (double)(fragNum * fragLength) / refLength, K_MER);
    printf("에러율 설정: %d%% | 실제 주입된 미스매치: %d개 / %ld bp (%.2f%%)\n\n",
           ERROR_RATE_PERCENT, injectedErrors, totalBases,
           totalBases > 0 ? (100.0 * injectedErrors / totalBases) : 0.0);

    // 공통 전처리: 인덱스 + 힙
    CountingIndex* countingIndex = buildCountingIndex(frags);
    MaxHeap* seedHeap = (countingIndex != NULL) ? buildSeedHeap(countingIndex) : NULL;

    if (countingIndex == NULL || seedHeap == NULL) {
        printf("인덱스/힙 생성에 실패했습니다.\n");
        freeMaxHeap(seedHeap);
        freeCountingIndex(countingIndex);
        for (int i = 0; i < fragNum; i++) free(frags[i]);
        free(frags);
        free(ref);
        return 1;
    }

    printTopSeeds(seedHeap, 5);

    // ===== [방식 1] 기존 greedy 조립 (벤치마크) =====
    clock_t s1 = clock();
    char* greedy = assembleReads(countingIndex, seedHeap, frags, MAX_MISMATCH);
    double t1 = (double)(clock() - s1) / CLOCKS_PER_SEC;

    // ===== [방식 2] Consensus 보정 조립 (우리 개선안) =====
    clock_t s2 = clock();
    char* consensus = assembleConsensus(countingIndex, seedHeap, frags, MAX_MISMATCH);
    double t2 = (double)(clock() - s2) / CLOCKS_PER_SEC;

    if (greedy != NULL) {
        printf("########## [방식 1] 기존 Greedy 조립 ##########");
        printPerformanceReport(ref, greedy, t1, calcMemory(countingIndex, seedHeap, strlen(greedy)));
    }
    if (consensus != NULL) {
        printf("\n########## [방식 2] Consensus(다수결 투표) 보정 조립 ##########");
        printPerformanceReport(ref, consensus, t2, calcMemory(countingIndex, seedHeap, strlen(consensus)));
    }

    // 메모리 해제
    free(consensus);
    free(greedy);
    freeMaxHeap(seedHeap);
    freeCountingIndex(countingIndex);
    for (int i = 0; i < fragNum; i++) free(frags[i]);
    free(frags);
    free(ref);

    return 0;
}
