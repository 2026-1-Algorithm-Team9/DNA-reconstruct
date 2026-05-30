#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "preIndex_core.h"

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

char** makeFrag(char* ref) {
    char** frags = (char**)malloc(fragNum * sizeof(char*));
    if (frags == NULL) return NULL;

    for (int i = 0; i < fragNum; i++) {
        frags[i] = (char*)malloc((fragLength + 1) * sizeof(char));

        int startIndex = rand() % (refLength - fragLength + 1);   // 원본에서 랜덤 시작 위치
        for (int j = 0; j < fragLength; j++) {
            frags[i][j] = ref[startIndex + j];                    // 원본의 substring을 그대로 떼옴
        }
        frags[i][fragLength] = '\0';
    }

    return frags;
}

int main(void) {
    char* ref = makeRef();
    char** frags = (ref != NULL) ? makeFrag(ref) : NULL;

    if (ref == NULL || frags == NULL) {
        printf("초기 데이터 생성에 실패했습니다.\n");
        free(frags);
        free(ref);
        return 1;
    }

    printf("=== 데이터 생성 파라미터 ===\n");
    printf("원본 길이(N): %d | 리드 길이(L): %d | 리드 개수(M): %d | 커버리지: %.1f배 | K-mer: %d\n\n",
           refLength, fragLength, fragNum,
           (double)(fragNum * fragLength) / refLength, K_MER);

    // ===== 알고리즘 시작 (시간 측정 구간: 데이터 생성/출력은 제외) =====
    clock_t start = clock();

    CountingIndex* countingIndex = buildCountingIndex(frags);
    MaxHeap* seedHeap = (countingIndex != NULL) ? buildSeedHeap(countingIndex) : NULL;

    char* assembled = NULL;
    if (countingIndex != NULL && seedHeap != NULL) {
        assembled = assembleReads(countingIndex, seedHeap, frags, MAX_MISMATCH);
    }

    clock_t end = clock();
    double duration = (double)(end - start) / CLOCKS_PER_SEC;
    // ===== 알고리즘 끝 =====

    if (countingIndex == NULL || seedHeap == NULL || assembled == NULL) {
        printf("알고리즘 수행에 실패했습니다.\n");
        free(assembled);
        freeMaxHeap(seedHeap);
        freeCountingIndex(countingIndex);
        for (int i = 0; i < fragNum; i++) free(frags[i]);
        free(frags);
        free(ref);
        return 1;
    }

    // 메모리 사용량 합산 (malloc 기준 간단 측정 = 공간복잡도 근사치)
    size_t memoryBytes = 0;
    memoryBytes += countingIndexMemory(countingIndex);                                 // 인덱스
    memoryBytes += maxHeapMemory(seedHeap);                                            // 힙
    memoryBytes += (size_t)fragNum * sizeof(char*) + (size_t)fragNum * (fragLength + 1); // 리드 조각들
    memoryBytes += (size_t)(refLength + 1);                                            // 원본 서열
    memoryBytes += (size_t)(fragNum * fragLength * 2 + 1);                             // 조립 임시 버퍼(peak)
    memoryBytes += strlen(assembled) + 1;                                              // 최종 조립 서열

    // 상위 시드 표시 (데이터가 작을 때만 상세 인덱스까지 출력)
    if (fragNum <= 12) {
        printCountingIndex(countingIndex, frags);
    }
    printTopSeeds(seedHeap, 5);

    // 5단계 결과 + 성능 리포트 (정확도/속도/메모리)
    printPerformanceReport(ref, assembled, duration, memoryBytes);

    // 메모리 해제
    free(assembled);
    freeMaxHeap(seedHeap);
    freeCountingIndex(countingIndex);
    for (int i = 0; i < fragNum; i++) free(frags[i]);
    free(frags);
    free(ref);

    return 0;
}
