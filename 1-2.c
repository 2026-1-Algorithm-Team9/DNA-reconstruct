#include <stdio.h>
#include <stdlib.h>

#define refLength 100   // 원본 염기서열 길이
#define fragLength 10   // 조각(Read) 길이
#define fragNum 10      // 조각 개수
#define K_MER 3         // 자를 k-mer의 길이
#define HASH_SIZE (1 << (K_MER * 2))    // 해시 크기

/* charToBit, freqArray만 있는 파일입니다 */

// 문자를 2비트 정수로 매핑 (A=00, C=01, G=10, T=11)
int charToBit(char c) {
    switch (c) {
        case 'A': return 0; 
        case 'C': return 1; 
        case 'G': return 2; 
        case 'T': return 3;
    }
}

int* freqArray(char** frags) {
    // 빈도 카운팅 배열 동적 할당 및 0으로 깨끗하게 초기화 (calloc 사용)
    int* countArray = (int*)calloc(HASH_SIZE, sizeof(int));
    if (countArray == NULL) return NULL;

    // k-mer 크기만큼만 비트를 남기기 위한 변수 cut
    int cut = (1 << (K_MER * 2)) - 1;

    for (int i = 0; i < fragNum; i++) {
        char* read = frags[i];
        int currentHash = 0;

        // 최초 k-mer 해시 생성 (첫 k글자 비트 조립)
        for (int j = 0; j < K_MER; j++) {
            currentHash = (currentHash << 2) | charToBit(read[j]);
        }
        // 첫 번째 해시 주소 방에 카운트 +1
        countArray[currentHash]++; 

        // 라빈-카프 기반 Rolling Hash로 슬라이딩
        for (int j = K_MER; j < fragLength; j++) {
            // 자릿수 밀고, 새 글자 합치고, cut으로 자르기
            currentHash = ((currentHash << 2) | charToBit(read[j])) & cut;
            
            // 튀어나온 해시 주소 방에 실시간 카운트 +1
            countArray[currentHash]++; 
        }
    }

    return countArray; // 계수 정렬 위한 완성된 빈도 수첩 반환
}