#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

// ===== 데이터 생성 파라미터 (조립 코드 preIndex_core.h와 통일) =====
#define refLength   20000                                  // 원본 게놈 길이 (N)
#define COVERAGE    60                                     // 목표 커버리지 (배수)
#define fragLength  100                                    // 리드 길이 (L)
#define fragNum     ((refLength * COVERAGE) / fragLength)  // 리드 개수 (M)
#define ERROR_RATE_PERCENT 1                               // 시퀀싱 에러율(%) — 0이면 에러 없음
#define MISMATCH_THRESHOLD 2                               // 매핑 시 허용 미스매치 수
#define K_MER_SIZE 4                                       // 시드로 쓸 k-mer 길이

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

static char mutateBase(char original) {
    char basis[4] = { 'A', 'C', 'G', 'T' };
    char c;
    do { c = basis[rand() % 4]; } while (c == original);
    return c;
}

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

// 문자를 0~3 숫자로 변환 (A,C,G,T 2비트 인코딩)
int char_to_int(char c) {
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'T') return 3;
    return 0;
}

// k-mer 해시값 (Rolling Hash, 비트 연산)
int get_hash(const char* str, int k) {
    int hash = 0;
    for (int i = 0; i < k; i++) hash = (hash << 2) | char_to_int(str[i]);
    return hash;
}

// ──────────────────────────────────────────────────────────────
// 해시 인덱스 매핑: 게놈을 k-mer 인덱스(Counting Sort)로 만들어,
// read의 첫 k-mer를 가진 후보 위치에서만 미스매치 검사 → 완전탐색보다 빠름.
// 반환: read가 매칭된 위치(미스매치 최소). 없으면 -1.
// ※ 게놈 인덱스는 read마다 새로 만들지 않고 main에서 한 번만 만들어 넘긴다.
// ──────────────────────────────────────────────────────────────
typedef struct {
    int* count;        // k-mer 빈도
    int* prefix_sum;   // 누적합(시작 좌표)
    int* index_array;  // k-mer별 게놈 위치 목록
    int hash_size;
} GenomeIndex;

GenomeIndex* build_index(const char* genome, int k) {
    int g_len = strlen(genome);
    int hash_size = 1 << (2 * k);

    GenomeIndex* gi = (GenomeIndex*)malloc(sizeof(GenomeIndex));
    gi->hash_size = hash_size;
    gi->count = (int*)calloc(hash_size + 1, sizeof(int));
    gi->prefix_sum = (int*)calloc(hash_size + 1, sizeof(int));

    for (int i = 0; i <= g_len - k; i++) gi->count[get_hash(genome + i, k)]++;
    for (int i = 1; i <= hash_size; i++)
        gi->prefix_sum[i] = gi->prefix_sum[i - 1] + gi->count[i - 1];

    int total_kmers = g_len - k + 1;
    gi->index_array = (int*)malloc(total_kmers * sizeof(int));
    int* cur = (int*)calloc(hash_size, sizeof(int));
    for (int i = 0; i <= g_len - k; i++) {
        int h = get_hash(genome + i, k);
        gi->index_array[gi->prefix_sum[h] + cur[h]] = i;
        cur[h]++;
    }
    free(cur);
    return gi;
}

void free_index(GenomeIndex* gi) {
    if (!gi) return;
    free(gi->count); free(gi->prefix_sum); free(gi->index_array); free(gi);
}

int hashing_search(const char* genome, const GenomeIndex* gi,
                   const char* read, int k, int threshold) {
    int g_len = strlen(genome);
    int r_len = strlen(read);
    if (k > r_len) return -1;

    int read_hash = get_hash(read, k);
    int start = gi->prefix_sum[read_hash];
    int count_match = gi->count[read_hash];

    int bestPos = -1, bestMismatch = threshold + 1;
    for (int i = 0; i < count_match; i++) {
        int g_idx = gi->index_array[start + i];
        if (g_idx + r_len > g_len) continue;
        int mismatch = 0;
        for (int j = 0; j < r_len; j++) {
            if (genome[g_idx + j] != read[j]) {
                mismatch++;
                if (mismatch >= bestMismatch) break;
            }
        }
        if (mismatch < bestMismatch) { bestMismatch = mismatch; bestPos = g_idx; }
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

    printf("========== Simple Hash 인덱스 매핑 (벤치마크) ==========\n");
    printf("[데이터] N=%d, L=%d, M=%d, 커버리지=%.1f배, 에러율=%d%%, k-mer=%d\n",
           refLength, fragLength, fragNum,
           (double)(fragNum * fragLength) / refLength, ERROR_RATE_PERCENT, K_MER_SIZE);

    // ===== 측정: 인덱스 구축 + 탐색 시간 + 매핑 정확도 =====
    int correct = 0, found = 0;
    clock_t s = clock();
    GenomeIndex* gi = build_index(genome, K_MER_SIZE);   // 게놈 인덱스 1회 구축
    for (int i = 0; i < fragNum; i++) {
        int pos = hashing_search(genome, gi, reads[i], K_MER_SIZE, MISMATCH_THRESHOLD);
        if (pos != -1) {
            found++;
            if (pos == trueStart[i]) correct++;
        }
    }
    double duration = (double)(clock() - s) / CLOCKS_PER_SEC;
    free_index(gi);

    double mapRate = 100.0 * found / fragNum;
    double accRate = 100.0 * correct / fragNum;

    printf("\n================ [Simple Hash] 성능 분석 리포트 ================\n");
    printf("[정확도] 매핑 성공률 : %6.2f %%   (%d / %d 리드)\n", mapRate, found, fragNum);
    printf("[정확도] 위치 정확도 : %6.2f %%   (정답 위치 %d / %d)\n", accRate, correct, fragNum);
    printf("[속도]   탐색 시간   : %.6f 초  (인덱스 구축 포함)\n", duration);
    printf("==============================================================\n");

    for (int i = 0; i < fragNum; i++) free(reads[i]);
    free(reads); free(genome); free(trueStart);
    return 0;
}
