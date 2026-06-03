#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

// ===== 데이터 생성 파라미터 (조립 코드 preIndex_core.h와 통일) =====
#define refLength   20000                                  // 원본 게놈 길이 (N)
#define COVERAGE    60                                     // 목표 커버리지 (배수)
#define fragLength  100                                    // 리드 길이 (L)
#define fragNum     ((refLength * COVERAGE) / fragLength)  // 리드 개수 (M)
#define ERROR_RATE_PERCENT 1                               // 시퀀싱 에러율(%) — 0이면 에러 없음
#define MISMATCH_THRESHOLD 2                               // 매핑 시 허용 미스매치 수

// OS와 무관하게 충분히 큰 난수 (Windows의 작은 RAND_MAX 보정)
static long bigRand(void) {
    return ((long)rand() << 15) | (long)rand();
}

char* makeRef() {
    char basis[4] = { 'A', 'T', 'C', 'G' };
    char* ref = (char*)malloc((refLength + 1) * sizeof(char));
    if (ref == NULL) return NULL;
    srand((unsigned int)time(NULL));
    for (int i = 0; i < refLength; i++) ref[i] = basis[rand() % 4];
    ref[refLength] = '\0';
    return ref;
}

// 한 염기를 자신과 다른 3개 중 하나로 치환 (미스매치 1개)
static char mutateBase(char original) {
    char basis[4] = { 'A', 'C', 'G', 'T' };
    char c;
    do { c = basis[rand() % 4]; } while (c == original);
    return c;
}

// 원본에서 substring을 떼와 리드 생성 + 에러 주입.
// trueStart[i] = i번째 리드가 원본에서 떼어진 정답 위치(매핑 정확도 채점용).
char** makeFrag(char* ref, int* trueStart) {
    char** frags = (char**)malloc(fragNum * sizeof(char*));
    if (frags == NULL) return NULL;
    for (int i = 0; i < fragNum; i++) {
        frags[i] = (char*)malloc((fragLength + 1) * sizeof(char));
        int startIndex = bigRand() % (refLength - fragLength + 1);
        if (trueStart) trueStart[i] = startIndex;
        for (int j = 0; j < fragLength; j++) {
            char base = ref[startIndex + j];
            if (ERROR_RATE_PERCENT > 0 && (rand() % 100) < ERROR_RATE_PERCENT)
                base = mutateBase(base);
            frags[i][j] = base;
        }
        frags[i][fragLength] = '\0';
    }
    return frags;
}

// Trivial 매핑 (완전 탐색): 게놈의 모든 위치를 한 칸씩 밀며 1:1 대조.
// 반환: read가 매칭된 위치(가장 미스매치가 적은 위치). 없으면 -1.
int trivial_search(const char* genome, const char* read, int threshold) {
    int g_len = strlen(genome);
    int r_len = strlen(read);
    int bestPos = -1, bestMismatch = threshold + 1;

    for (int i = 0; i <= g_len - r_len; i++) {
        int mismatch = 0;
        for (int j = 0; j < r_len; j++) {
            if (genome[i + j] != read[j]) {
                mismatch++;
                if (mismatch >= bestMismatch) break;   // 현재 최선보다 나쁘면 중단
            }
        }
        if (mismatch < bestMismatch) { bestMismatch = mismatch; bestPos = i; }
        if (bestMismatch == 0) break;                  // 완벽 매칭이면 더 볼 것 없음
    }
    return (bestMismatch <= threshold) ? bestPos : -1;
}

int main(void) {
    int* trueStart = (int*)malloc(fragNum * sizeof(int));
    char* genome = makeRef();
    char** reads = (genome != NULL) ? makeFrag(genome, trueStart) : NULL;
    if (genome == NULL || reads == NULL || trueStart == NULL) {
        printf("데이터 생성 실패\n"); return 1;
    }

    printf("========== Trivial 완전탐색 매핑 (벤치마크) ==========\n");
    printf("[데이터] N=%d, L=%d, M=%d, 커버리지=%.1f배, 에러율=%d%%\n",
           refLength, fragLength, fragNum,
           (double)(fragNum * fragLength) / refLength, ERROR_RATE_PERCENT);

    // ===== 측정: 탐색 시간 + 매핑 정확도 =====
    int correct = 0, found = 0;
    clock_t s = clock();
    for (int i = 0; i < fragNum; i++) {
        int pos = trivial_search(genome, reads[i], MISMATCH_THRESHOLD);
        if (pos != -1) {
            found++;
            if (pos == trueStart[i]) correct++;   // 정답 위치를 맞췄나
        }
    }
    double duration = (double)(clock() - s) / CLOCKS_PER_SEC;

    double mapRate = 100.0 * found / fragNum;       // 매칭을 찾아낸 비율
    double accRate = 100.0 * correct / fragNum;     // 정답 위치까지 맞춘 비율

    printf("\n================ [Trivial] 성능 분석 리포트 ================\n");
    printf("[정확도] 매핑 성공률 : %6.2f %%   (%d / %d 리드)\n", mapRate, found, fragNum);
    printf("[정확도] 위치 정확도 : %6.2f %%   (정답 위치 %d / %d)\n", accRate, correct, fragNum);
    printf("[속도]   탐색 시간   : %.6f 초\n", duration);
    printf("==========================================================\n");

    for (int i = 0; i < fragNum; i++) free(reads[i]);
    free(reads); free(genome); free(trueStart);
    return 0;
}
