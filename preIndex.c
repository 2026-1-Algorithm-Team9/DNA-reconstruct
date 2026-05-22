#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define refLength 100   // 원본 염기서열 길이
#define fragLength 10   // 조각(Read) 길이
#define fragNum 10      // 조각 개수
#define K_MER 3         // 자를 k-mer의 길이
#define HASH_SIZE (1 << (K_MER * 2))    // 해시 크기

/* charToBit, freqArray 함수 구현했습니다
   charToBit, freqArray 함수랑 main함수 일부만 가져다 쓰면 될 것 같아용
   main에는 charToBit, freqArray 함수가 잘 동작하는지 확인용 코드를 적어놓아서 나중에 지울 것들 지우면 될듯요 */

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

    int* freqResult = freqArray(frags);

    // 해시의 몇 번 방에 카운트가 쌓였는지 확인 코드
    printf("=== 빈도 카운팅 정적 배열 (%d개의 방) ===\n", HASH_SIZE);
    for (int i = 0; i < HASH_SIZE; i++) {
        if (freqResult[i] > 0) {
            printf("해시 주소 [%2d번 방] -> 등장 빈도: %d회\n", i, freqResult[i]);
        }
    }

    // 메모리 해제
    free(freqResult);
    for (int i = 0; i < fragNum; i++) free(frags[i]);
    free(frags);
    free(ref);

    return 0;
}