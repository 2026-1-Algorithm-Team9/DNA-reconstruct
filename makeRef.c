#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define refLength 100   // 원본 염기서열 길이
#define fragLength 10   // 조각 길이
#define fragNum 10      // 조각 개수

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
    char* ref = makeRef();        // refLength개의 염기서열을 생성하여 ref에 저장 -> 염기 서열의 길이 변경 가능

    printf("<Reference Sequence>\n");
    printf("%s\n\n", ref);

    char** frags = makeFrag(ref);     // ref에서 fragNum개의 조각을 fragLength 길이로 만들어서 frags에 저장

    for (int i = 0; i < fragNum; i++) {     // 조각들이 잘 저장되었는지 확인용.. 나중에 지워도 됨
        printf("%s\n", frags[i]);
    }

    for (int i = 0; i < fragNum; i++) {
        free(frags[i]); // 개별 조각 해제
    }
    free(frags); // 조각 포인터 배열 해제
    free(ref);   // 원본 서열 해제

    return 0;
}