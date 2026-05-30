#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "preIndex_core.h"

#define refLength 1000   // 원본 염기서열 길이

char* makeRef() {
    char basis[4] = {'A', 'T', 'C', 'G'};

    // 염기 서열 저장하기 위한 배열 동적 할당 (+1은 문자열 종료 문자 '\0'를 위해)
    char* ref = (char*)malloc((refLength + 1) * sizeof(char));

    for (int i = 0; i < refLength; i++) {
        int index = rand() % 4;
        ref[i] = basis[index];
    }
    ref[refLength] = '\0'; // 문자열 종료 문자 추가

    return ref;
}

// preIndex.c 의 makeFrag 함수
char** makeFrag(char* ref) { // 원본 서열을 매개변수로 받아서
    char** frags = (char**)malloc(fragNum * sizeof(char*)); 

    for (int i = 0; i < fragNum; i++) {
        frags[i] = (char*)malloc((fragLength + 1) * sizeof(char)); 

        // 원본 중 어디서 자를지 '시작 위치'를 랜덤으로 결정
        int startIndex = rand() % (refLength - fragLength + 1); 

        for (int j = 0; j < fragLength; j++) {
            // 글자 내용은 원본에서 떼어오지만, 순서는 완전히 뒤죽박죽 섞이게 됩니다.
            frags[i][j] = ref[startIndex + j]; 
        }
        frags[i][fragLength] = '\0'; 
    }
    return frags;
}

int main(void) {
    clock_t total_start = clock();
    srand(time(NULL));

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


    // ===============================================================
    // 조립 및 출력 구간
    // ===============================================================
    if (countingIndex != NULL && seedHeap != NULL) {
        // 조립 함수 실행
        char* assembledRef = assembleReads(countingIndex, seedHeap, frags, MAX_MISMATCH);
        
        // 메모리 해제 직전, 콘솔 출력이 모두 끝난 시점에서 종료 시간 기록
        clock_t total_end = clock();
        double total_duration = (double)(total_end - total_start) / CLOCKS_PER_SEC;

        // 최종 결과물 예쁘게 출력하기 (전체 소요 시간을 넘겨줍니다)
        printPerformanceReport(ref, assembledRef, total_duration);

        // 메모리 해제
        if (assembledRef != NULL) free(assembledRef);
    }


    // 메모리 해제
    freeMaxHeap(seedHeap);
    freeCountingIndex(countingIndex);
    for (int i = 0; i < fragNum; i++) free(frags[i]);
    free(frags);
    free(ref);

    return 0;
}