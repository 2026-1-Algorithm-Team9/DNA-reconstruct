#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

#define refLength 100   // 원본 염기서열 길이
#define fragLength 10   // 조각 길이
#define fragNum 10      // 조각 개수

char* makeRef() {
    char basis[4] = { 'A', 'T', 'C', 'G' };

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
        frags[i] = (char*)malloc((fragLength + 1) * sizeof(char));     // 각 조각을 저장하기 위한 배열 동적 할당 (+1은 문자열 종료 문자 '\0'를 위해)

        int startIndex = rand() % (refLength - fragLength + 1);     // 원본 염기서열에서 조각을 시작할 랜덤한 인덱스 생성
        for (int j = 0; j < fragLength; j++) {
            frags[i][j] = ref[startIndex + j];     // 원본 염기서열에서 조각을 추출하여 frags 배열에 저장
        }
        frags[i][fragLength] = '\0';      // 문자열 종료 문자 추가
    }

    return frags;
}

// Trivial 매핑 알고리즘 (완전 탐색)
// genome: 원본 참조 서열, read: 찾고자 하는 짧은 서열, threshold: 허용되는 최대 미스매치 수

void trivial_search(const char* genome, const char* read, int threshold) {

    int g_len = strlen(genome);
    int r_len = strlen(read);

    // 원본 서열의 처음부터 끝까지 한 칸씩 이동하며 대조
    for (int i = 0; i <= g_len - r_len; i++) {

        int mismatch = 0;

        // read 길이만큼 문자열을 1:1로 비교
        for (int j = 0; j < r_len; j++) {

            if (genome[i + j] != read[j]) {

                mismatch++;

                // 허용치(threshold)를 넘어가면 즉시 비교 중단 (최적화)
                if (mismatch > threshold) break;

            }

        }

        // 미스매치 수가 허용치 이내라면 결과 출력
        if (mismatch <= threshold) {

            printf("Match at index %d (mismatches: %d)\n", i, mismatch);

        }

    }

}

int main() {

    // 테스트용 시퀀스 설정
    char* genome = makeRef();
    char** reads = makeFrag(genome);

    int mismatch_threshold = 1; // 1개까지 불일치 허용

    printf("원본 서열: %s\n\n", genome);

    // 10개의 조각에 대해 탐색 실행
    for (int i = 0; i < fragNum; i++) {
        printf("--- [%d번째 조각: %s] Trivial Search ---\n", i + 1, reads[i]);
        trivial_search(genome, reads[i], mismatch_threshold);
        printf("\n");
    }

    // 메모리 해제
    for (int i = 0; i < fragNum; i++) {
        free(reads[i]); // 개별 조각 해제
    }
    free(reads);  // 조각 포인터 배열 해제
    free(genome); // 원본 서열 해제

    return 0;

}