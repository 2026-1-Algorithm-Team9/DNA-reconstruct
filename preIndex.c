#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "preIndex_core.h"

#define refLength 100   // 원본 염기서열 길이

char* makeRef() {
    char basis[4] = {'A', 'T', 'C', 'G'};

    // 염기 서열 저장하기 위한 배열 동적 할당 (+1은 문자열 종료 문자 '\0'를 위해)
    char* ref = (char*)malloc((refLength + 1) * sizeof(char));
    srand(time(NULL));

    for (int i = 0; i < refLength; i++) {
        int index = rand() % 4;
        ref[i] = basis[index];
    }
    ref[refLength] = '\0'; // 문자열 종료 문자 추가

    return ref;
}

char** makeFrag(char* ref) {
    char** frags = (char**)malloc(fragNum * sizeof(char*));     // 조각들을 저장하기 위한 2차원 배열을 동적으로 할당

    for (int i = 0; i < fragNum; i++) {
        frags[i] = (char*)malloc((fragLength+1) * sizeof(char));     // 각 조각을 저장하기 위한 배열 동적 할당 (+1은 문자열 종료 문자 '\0'를 위해)

        int startIndex = rand() % (refLength - fragLength + 1);     // 원본 염기서열에서 조각을 시작할 랜덤한 인덱스 생성
        for (int j = 0; j < fragLength; j++) {
            frags[i][j] = ref[startIndex + j];     // 원본 염기서열에서 조각을 추출하여 frags 배열에 저장
        }
        frags[i][fragLength] = '\0';      // 문자열 종료 문자 추가
    }

    return frags;
}

int main(void) {
    char* ref = makeRef();
    char** frags = makeFrag(ref);  

    printf("=== 생성된 원본 염기 서열 ===\n");
    printf("%s\n\n", ref);

    printf("=== 생성된 숏 리드 조각 %d개 ===\n", fragNum);
    for (int i = 0; i < fragNum; i++) printf("[%d] %s\n", i + 1, frags[i]);
    printf("\n");

    CountingIndex* countingIndex = buildCountingIndex(frags);
    MaxHeap* seedHeap = NULL;

    if (countingIndex != NULL) {
        seedHeap = buildSeedHeap(countingIndex);
    }

    if (countingIndex == NULL || seedHeap == NULL) {
        printf("메모리 할당에 실패했습니다.\n");

        freeMaxHeap(seedHeap);
        freeCountingIndex(countingIndex);
        for (int i = 0; i < fragNum; i++) free(frags[i]);
        free(frags);
        free(ref);
        return 1;
    }

    // 해시의 몇 번 방에 카운트가 쌓였는지 확인 코드
    printf("=== 빈도 카운팅 정적 배열 (%d개의 방) ===\n", HASH_SIZE);
    for (int i = 0; i < HASH_SIZE; i++) {
        if (countingIndex->countArray[i] > 0) {
            printf("해시 주소 [%2d번 방] -> 등장 빈도: %d회\n", i, countingIndex->countArray[i]);
        }
    }
    printf("\n");

    printCountingIndex(countingIndex, frags);
    printTopSeeds(seedHeap, 5);

    // 메모리 해제
    freeMaxHeap(seedHeap);
    freeCountingIndex(countingIndex);
    for (int i = 0; i < fragNum; i++) free(frags[i]);
    free(frags);
    free(ref);

    return 0;
}